/***************************************************************************

Copyright (c) 2025 Kay F. Jahnke

Licensed under MIT License:

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

****************************************************************************/

// this header defines the spec_t class - this class holds the
// parameterization for a rendering job - and gives access to the
// rendered data. spec_t objects are created by 'visor' and pushed
// to a deque in shared memory. 'render jobs' pull them from the deque,
// render the frames and then push the job to the frame queue, now
// 'loaded' with frame data.
// This header also codes access and structure of the shared memory
// which is used to run the 'visor protocol'. The protocol itself is
// partly coded in render.cc, and partly in static member functions
// of class ipc_data_t, to make using the visor protocol simple for
// rendering processes: all they need to do is to #include visor.h
// and then pass a job handling functor to the render_loop member
// function. The design expects that the visor process is used with
// varying render processes - further commodification on the visor
// visor side would be feasible as well.

#include <cstddef>
#include <string>
#include <iostream>
#include <chrono>
#include <thread>
#include <functional>

#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/containers/deque.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/containers/string.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>
#include <boost/interprocess/sync/interprocess_mutex.hpp>
#include <boost/interprocess/sync/named_mutex.hpp>
#include <boost/interprocess/sync/interprocess_condition.hpp>
#include <boost/interprocess/offset_ptr.hpp>

typedef std::chrono::high_resolution_clock::time_point vs_time_t ;
#define NOW() std::chrono::high_resolution_clock::now()
#define delta_t(a,b) std::chrono::duration_cast<std::chrono::milliseconds> \
                       ( b - a ) . count() ;

const std::size_t NJOBS = 12 ;
const std::size_t NFRAMES = 5 ;

using namespace boost::interprocess;

using int_allocator = managed_shared_memory::allocator<int>::type;

// Typ-Alias für den Offset-Pointer, der auf rohen Speicher zeigt
using ShmemBytePtr = boost::interprocess::offset_ptr<std::byte>;

struct spec_t
{
  std::size_t serial_no ;
  bool refine = false ;
  bool snapshot = false ;

  boost::interprocess::string filename ;

  int prj_img ;
  double hfov_img ;
  double vfov_img ;
  double yaw_img ;
  double pitch_img ;
  double roll_img ;

  int prj_cam ;
  double hfov_cam ;
  double vfov_cam ;
  double yaw_cam ;
  double pitch_cam ;
  double roll_cam ;
  unsigned int width_cam ;
  unsigned int height_cam ;
  float brighten ;

  int bits_per_pixel = 32 ;
  std::size_t buffer_index ;

  vs_time_t job_enqueued ;
  vs_time_t rendering_starts ;
  vs_time_t rendering_ends ;
  vs_time_t frame_enqueued ;
  vs_time_t frame_fetched ;
  vs_time_t texture_updated ;
  vs_time_t sprite_created ;
  vs_time_t display_called ;
  vs_time_t display_returns ;

  void print_timing()
  {
    auto dt1 = delta_t ( job_enqueued , rendering_starts ) ;
    auto dt2 = delta_t ( rendering_starts , rendering_ends ) ;
    auto dt3 = delta_t ( rendering_ends , frame_enqueued ) ;
    auto dt4 = delta_t ( frame_enqueued , frame_fetched ) ;
    auto dt5 = delta_t ( frame_fetched , texture_updated ) ;
    auto dt6 = delta_t ( texture_updated , sprite_created ) ;
    auto dt6a = delta_t ( sprite_created , display_called ) ;
    auto dt7 = delta_t ( display_called , display_returns ) ;
    auto dt8 = delta_t ( job_enqueued , display_returns ) ;

    std::cout << "--------------" << serial_no << std::endl ;
    std::cout << "rnd_starts:   " << dt1 << std::endl ;
    std::cout << "rendering:    " << dt2 << std::endl ;
    std::cout << "to frame q:   " << dt3 << std::endl ;
    std::cout << "from frame q: " << dt4 << std::endl ;
    std::cout << "to texture:   " << dt5 << std::endl ;
    std::cout << "to sprite:    " << dt6a << std::endl ;
    std::cout << "to call dsp:  " << dt6 << std::endl ;
    std::cout << "to dsp ret:   " << dt7 << std::endl ;
    std::cout << "total:        " << dt8 << std::endl ;
  }
} ;

using spec_t_allocator = managed_shared_memory::allocator<spec_t>::type;
using spec_t_array = boost::interprocess::vector<spec_t, spec_t_allocator>;

struct info_t
{
  std::size_t desktop_width ;
  std::size_t desktop_height ;
  std::size_t nframes ;
  ShmemBytePtr data_ofsptr = nullptr ;
  ShmemBytePtr get_buffer_address ( std::size_t i )
  {
    if ( i < NFRAMES )
      return data_ofsptr + desktop_width * desktop_height * 4U * i ;
    return nullptr ;
  }
  void init()
  {
    for ( std::size_t i = 0 ; i < NFRAMES ; i++ )
    {
      auto * p = (int*) get_buffer_address ( i ) . get() ;
      for ( std::size_t j = 0 ; j < desktop_width * desktop_height ; j++ )
        p[j] = 0xFFAA0000 | ( i * 50 + j ) ;
    }
  }
} ;

// store_t implements a simple stack of NFRAMES std::size_t.
// the get and put member functions return true on success, on
// failure they return false and there is no effect. The stack
// used to handle inidices pertaining to the set of frame buffers
// held in shared memory. These buffers are sized to contain as
// many pixels as a full desktop, the actual content may be less
// if the buffers are used to contain texture data for smaller
// frames, but there is no way to store larger frames. This is
// a compromise: we don't want to allocate new frame buffers if
// the window size changes, so we need a fixed limit, and we
// assume that the mechanism is used to provide screen data.

struct store_t
{
  boost::interprocess::interprocess_mutex mutex ;
  std::size_t store [ NFRAMES ] ;
  std::size_t lid ;
  const std::size_t capacity ;

  store_t()
  : capacity ( NFRAMES )
  {
    for ( std::size_t i = 0 ; i < NFRAMES ; i++ )
      store[i] = NFRAMES - 1U - i ;
    lid = NFRAMES ;
  }

  bool get ( size_t & storee )
  {
    boost::interprocess::scoped_lock
      <boost::interprocess::interprocess_mutex> lock ( mutex ) ;

    if ( lid )
    {
      storee = store [ --lid ] ;
      return true ;
    }
    return false ;
  }

  bool put ( const size_t & storee )
  {
    boost::interprocess::scoped_lock
      <boost::interprocess::interprocess_mutex> lock ( mutex ) ;

    if ( lid < capacity )
    {
      store [ lid++ ] = storee ;
      return true ;
    }
    return false ;
  }
} ;

// for now, a very simple construct to serialize and de-serialize
// an argument vector. we use a fixed-size buffer inside the ipc_data
// structure and concatenate all arguments, including their trailing
// zero. We add another zero after the last arg.
// extraction does not actually extract the arguments but only pointers
// to the location where the arguments lie inside the buffer, which has
// very little overhead.

const std::size_t FLAT_ARGS_SZ = 8192 ;

struct flat_args_t
{
  char buffer [ FLAT_ARGS_SZ ] ;
  std::size_t wpos ;

  flat_args_t()
  : wpos ( 0U )
  { }

  bool append ( std::size_t argc , const char ** argv )
  {
    auto _wpos = wpos ;
    for ( std::size_t i = 0 ; i < argc ; i++ )
    {
      const char * p_src = argv[i] ;
      assert ( p_src != nullptr ) ;
      if ( *p_src )
      {
        while ( *p_src )
        {
          if ( wpos >= FLAT_ARGS_SZ - 3U )
          {
            // bail out returning wpos to safe location
            wpos = _wpos ;
            return false ;
          }
          buffer [ wpos++ ] = * p_src++ ;
        }
        buffer [ wpos++ ] = '\0' ;
      }
    }
    if ( _wpos != wpos )
    {
      buffer [ wpos++ ] = '\0' ;
    }
    return true ;
  }

  bool extract ( std::size_t & argc , std::vector < const char * > & argv )
  {
    argc = 0 ;
    argv.clear() ;
    if ( wpos < 2L )
      return false;
    assert ( buffer [ wpos - 1U ] == '\0' && buffer [ wpos - 2U ] == '\0' ) ;

    const char * p_src = buffer  ;
    while ( true )
    {
      argv.push_back ( p_src ) ;
      ++argc ;
      while ( *p_src++ ) { }
      if ( *p_src == '\0' )
        break ;
    }
    return true ;
  }
} ;

using segment_manager_t
 = boost::interprocess::managed_shared_memory::segment_manager;
using job_t_deque = boost::interprocess::deque<int, int_allocator> ;
using ipc_mutex_t = boost::interprocess::interprocess_mutex ;
using ipc_condition_t = boost::interprocess::interprocess_condition ;
using ipc_scoped_lock_t = boost::interprocess::scoped_lock<ipc_mutex_t> ;

struct ipc_data_t
{
  // first, some members for shared persistent parameters
  // - this was in struct info_t

  // extent of the desktop - we need this size to shape the chunks
  // of frame buffer memory

  std::size_t desktop_width ;
  std::size_t desktop_height ;

  // this pointer points to the first frame buffer

  ShmemBytePtr data_ofsptr = nullptr ;

  // obtain a pointer to frame buffer #i

  ShmemBytePtr get_buffer_address ( std::size_t i )
  {
    if ( i < NFRAMES )
      return data_ofsptr + desktop_width * desktop_height * 4U * i ;
    return nullptr ;
  }

  // optional - set the frame buffers to a defined initial state

  void init()
  {
    for ( std::size_t i = 0 ; i < NFRAMES ; i++ )
    {
      auto * p = (int*) get_buffer_address ( i ) . get() ;
      for ( std::size_t j = 0 ; j < desktop_width * desktop_height ; j++ )
        p[j] = 0xFFAA0000 | ( i * 50 + j ) ;
    }
  }

  // now some members which are used for IPC, like mutexes and condition
  // variables.

  job_t_deque job_queue ;
  job_t_deque frame_queue ;
  spec_t_array spec_array ;

  // ipc_data_t's c'tor takes a segment_manager_t* and passes it to the
  // data members which also need it in their c'tor

  ipc_data_t ( std::size_t _desktop_width ,
               std::size_t _desktop_height ,
               segment_manager_t * p_manager )
  : desktop_width ( _desktop_width ) ,
    desktop_height ( _desktop_height ) ,
    job_queue ( p_manager ) ,
    frame_queue ( p_manager ) ,
    spec_array ( NJOBS , p_manager )
  {
    std::size_t nbytes = NFRAMES * desktop_width * desktop_height * 4U ;
    data_ofsptr = static_cast<std::byte*> ( p_manager->allocate ( nbytes ) ) ;
    assert ( data_ofsptr != nullptr ) ;
  }

  // 'store' is a helper structure providing a stack of frame buffer inidices
  
  store_t store ;

  // flat_args contains a serialized argv

  flat_args_t flat_args ;

  // now some mutexes and cv's

  ipc_mutex_t job_mutex ;
  ipc_mutex_t frame_mutex ;
  ipc_mutex_t work_mutex ;

  ipc_condition_t job_cv ;
  ipc_condition_t frame_cv ;
  ipc_condition_t work_cv ;
} ;

struct visor_protocol
{
  // setting up and using the shared memory is coded here, in static
  // member functions. this is to make it easy to code processes which
  // use the visor protocol: shm_setup creates and populates the shared
  // memory used for IPC, and 'render_loop' handles all the interaction
  // with the shared memory to obtain jobs, and calls handle_job for
  // each job it obtains - handle_job is provided by the caller, which
  // can remain unaware of all the complexity.
  // note that pushing jobs to the job queue and obtaining frames is
  // coded in 'visor'. It would be feasible to also commodify that code
  // and put it into this header. For now we stick with a single visor
  // process and potentially several flavours of render process - e.g.
  // there is one 'dummy' render process which only serves as 'sparring
  // partner' and doesn't yield useful data.

  static bool shm_setup ( std::size_t max_frame_width ,
                          std::size_t max_frame_height )
  {
    try
    {
      // destroy visor_shm and the named mutex if it's present.

      boost::interprocess::shared_memory_object::remove ( "visor_shm" ) ;
      boost::interprocess::named_mutex::remove("visor_shm_ready_mutex");

      // create shared memory segment. first calculate the size:

      auto full_frame_size = 4U * max_frame_width * max_frame_height ;

      // we want five full frames to rotate

      auto need = NFRAMES * full_frame_size ;

      // and space for NJOBS spec_t objects

      need += NJOBS * sizeof ( spec_t ) ;

      // plus some space for mutexes, queues and condition variables
      // currently double 'need', due to duplicated infrastructure

      need += FLAT_ARGS_SZ + 20000U ;

      std::cout << "allocating " << need << " bytes in shm"
                << std::endl ;

      boost::interprocess::managed_shared_memory segment(
          boost::interprocess::create_only,
          "visor_shm", need ) ;

      // now we create a named mutex in the newly acquired shared memory
      // segment to prevent the rendering process from accessing it while
      // we're still setting up the ipc_data structure

      boost::interprocess::named_mutex shm_ready_mutex
        ( boost::interprocess::open_or_create , 
          "visor_shm_ready_mutex" ) ;

      // Lock the mutex

      boost::interprocess::scoped_lock<boost::interprocess::named_mutex> 
          lock(shm_ready_mutex);
      
      auto * p_manager = segment.get_segment_manager() ;

      ipc_data_t * p_ipc_data
        = segment.find_or_construct<ipc_data_t>("ipc_data")
            ( max_frame_width , max_frame_height , p_manager ) ;

      std::cout << "visor: constructed ipc object" << std::endl ;

      // the lock is released now that the scope closes, and the render
      // process can gain access to the shared memory
    }
    catch (const boost::interprocess::interprocess_exception& ex)
    {
      std::cerr << "error: " << ex.what() << std::endl;
      return false ;
    }
    return true ;
  }

  typedef std::function < bool ( ipc_data_t & , int ) > job_handler_t ;

  static bool render_loop ( job_handler_t job_handler )
  {
    // we code so that even if the render process is launched separately,
    // the connection via shared memory can still be established. So if
    // gaining access to the shared memory fails, we try again several
    // times until we give up. IF the visor process launches the render
    // process, gaining access to the shared memory will work immediately,
    // but coding for a possible delay gives us flexibility to handle
    // the conncetion either way.

    boost::interprocess::managed_shared_memory * p_shm ;

    using namespace std::chrono_literals ;
    auto retry_after = 20ms ;

    for ( int tries = 0 ; tries <= 10 ; tries++ )
    {
      if ( tries == 10 )
      {
        std::cerr << "could not gain access to visor's shared memory"
                  << std::endl ;
        exit ( -1 ) ;
      }

      if ( tries )
      {
        std::cout << "micro nap before trying again" << std::endl ;
        std::this_thread::sleep_for ( tries * retry_after ) ;
      }

      try
      {
        p_shm = new boost::interprocess::managed_shared_memory
          ( boost::interprocess::open_only , "visor_shm" ) ;

        break ;
      }
      catch (const boost::interprocess::interprocess_exception& ex)
      {
        std::cerr << "error: " << ex.what() << std::endl ;
      }
    }

    assert ( p_shm != nullptr ) ;

    std::cout << "render: have segment" << std::endl ;
    
    // we now try and obtain a lock on the named mutex which protects
    // te setting-up of the ipc_data structure. gaining the lock means
    // the data are ready - until then, the process blocks.
    {
      boost::interprocess::named_mutex shm_ready_mutex
        ( boost::interprocess::open_or_create , 
          "visor_shm_ready_mutex" ) ;
      
      // Lock the mutex

      boost::interprocess::scoped_lock<boost::interprocess::named_mutex> 
          lock(shm_ready_mutex);
    }

    std::pair<ipc_data_t*,std::size_t>
      ret = p_shm->find<ipc_data_t>("ipc_data");
    
    ipc_data_t & ipc ( * ( ret.first ) ) ;
    std::cout << "render: have ipc handle" << std::endl ;

    std::vector < const char * > visor_argv ;
    std::size_t visor_argc ;
    
    ipc.flat_args.extract ( visor_argc , visor_argv ) ;

    // for ( std::size_t i = 0 ; i < visor_argc ; i++ )
    //   std::cout << "visor arg " << i << " "
    //             << visor_argv [ i ] << std::endl ;

    int job_pending ;
    bool have_job ;

    while ( true )
    {
      have_job = false ;

      {
        // we protect the queue access with a 'scoped_lock' which is the
        // euqivalent of a std::lock_guard. if the job queue has a job
        // and the frame queue is not yet filled to the limit, we extract
        // the job by copying it to 'job_pending' and set a flag. The copy
        // is done to close the scope holding the lock as quickly as
        // possible, rather than handling the job while the lock is held.

        boost::interprocess::scoped_lock
          <boost::interprocess::interprocess_mutex> lock(ipc.job_mutex);

        if ( ipc.job_queue.size() )
        {
          job_pending = ipc.job_queue.front() ;
          ipc.job_queue.pop_front() ;
          have_job = true ;
        }
        else
        {
          // if there are no jobs in the queue, we wait on job_cv:
          // the main process will notify whenever it pushes a job
          // to the job queue.

          ipc.job_cv.wait ( lock ) ;
        }
      }

      if ( have_job )
      {
        // now we have a job. we handle it. what the queue yields is just
        // a job_t, which is an index into spec_array. we pass the spec_t
        // per reference.

        spec_t & spec ( ipc.spec_array [ job_pending ] ) ;

        // job with serial number zero signals session end

        if ( spec.serial_no == 0 )
          break ;

        spec.rendering_starts = NOW() ;
        bool snapshot = spec.snapshot ;
        bool success = job_handler ( ipc , job_pending ) ;
        spec.rendering_ends = NOW() ;

        if ( ! success )
          break ;

        // if the job was a snapshot, we mustn't push anything to the
        // frame queue - the frame was stored in a file and it's also
        // potentially much larger than one of the frame buffers in
        // shared memory.

        if ( snapshot )
          continue ;

        // we want to push a frame to the frame queue, but to do so, the
        // frame queue must have space. if the queue is full already, we
        // wait until the mainprocess signals that it has taken a frame
        // from the queue.

        while ( true )
        {
          {
            boost::interprocess::scoped_lock
              <boost::interprocess::interprocess_mutex> lock(ipc.frame_mutex);

            if ( ipc.frame_queue.size() < 3 )
            {
              spec.frame_enqueued = NOW() ;
              ipc.frame_queue.push_back ( job_pending ) ;
              have_job = false ;
              break ;
            }
            else
            {
              // the frame queue is already filled to the limit, we can't
              // yet push the frame we have made. we wait for the main
              // thread to notify on frame_cv, which it does whenever it
              // consumes a frame by pushing it to the GPU

              ipc.frame_cv.wait ( lock ) ;
            }
          }
        }

        // with or without a wait, we have succeeded in pushing the frame
        // to the frame queue. now we notify the main thread of the push,
        // in case it's waiting on work_cv.
        
        ipc.work_cv.notify_one() ;
      }
    }

    std::cout << "render job terminates" << std::endl ;

    // destroy the shared memory. this also destroys all data residing
    // inside the shared memory: queues, jobs, buffers. We needn't
    // bother with deallocating them.

    boost::interprocess::shared_memory_object::remove ( "visor_shm" ) ;

    // we also remove visor_shm_ready_mutex, which exists independently
    // of the managed shared memory object.

    boost::interprocess::named_mutex::remove("visor_shm_ready_mutex");

    delete p_shm ;
    return 0 ;
  }
} ;



