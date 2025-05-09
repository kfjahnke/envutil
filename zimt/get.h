/************************************************************************/
/*                                                                      */
/*    zimt - abstraction layer for SIMD programming                     */
/*                                                                      */
/*            Copyright 2023 by Kay F. Jahnke                           */
/*                                                                      */
/*    The git repository for this software is at                        */
/*                                                                      */
/*    https://github.com/kfjahnke/zimt                                  */
/*                                                                      */
/*    Please direct questions, bug reports, and contributions to        */
/*                                                                      */
/*    kfjahnke+zimt@gmail.com                                           */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

/*! \file get.h

    \brief classes to provide input to an 'act' functor

    zimt::process generates input to the 'act' functor, working along
    the aggregation axis, while the coordinate's component along the
    other axes remain constant. This effectively 1D operation can
    be coded efficiently, whereas the 'act' functor would typically
    look at it's input as an nD entity and not be aware of the fact
    that only one of the input's components is 'moving'.
    This header provides a set of classes which are meant to fit into
    the slot of generating input to the 'act' functor. They follow a
    specific pattern which results from the logic in zimt_process,
    namely the two init and two increase overloads.
*/

#if defined(ZIMT_GET_H) == defined(HWY_TARGET_TOGGLE)
  #ifdef ZIMT_GET_H
    #undef ZIMT_GET_H
  #else
    #define ZIMT_GET_H
  #endif

HWY_BEFORE_NAMESPACE() ;
BEGIN_ZIMT_SIMD_NAMESPACE(zimt)

// Here we have a collection of get_t objects to cover a set of
// common data acquisition strategies, and to serve as templates
// for your own creations. Note how all of these classes use the
// same set of template arguments (even if at times not all of
// them are actually used - this is for syntactic uniformity.
/// class get_crd is an implementation of a 'get_t' class which the
/// rolling-out code uses to produce input values to the functor 'act'.
/// This specific class provides discrete nD coordinates. The c'tor
/// receives the 'hot' axis along which the coordinate will vary,
/// the other components remain constant.

template < typename T ,    // elementary/fundamental type
           std::size_t N , // number of channels
           std::size_t D , // dimension of the view/array
           std::size_t L = vector_traits < T > :: vsize >
struct get_crd
{
  typedef zimt::xel_t < T , N > value_t ;
  typedef simdized_type < value_t , L > value_v ;
  typedef typename value_v::value_type value_ele_v ;
  typedef zimt::xel_t < long , D > crd_t ;

  const std::size_t d ;

  // get_crd's c'tor receives the processing axis

  get_crd ( const zimt::bill_t & bill = zimt::bill_t() )
  : d ( bill.axis )
  { }

  // init is used to initialize the vectorized value to the value
  // it should hold at the beginning of the run. The discrete
  // coordinate 'crd' gives the location of the first value, and
  // the get_crd infers the start value from it. Note that we can't
  // simply start with a value of zero for the 'hot' axis, because
  // zimt::process may 'cut up' the 1D subarrays further into
  // 'segments'. Since this is an archetypal implementation of a
  // get_t class, I'll explain the signature some more: The first
  // argument is a reference to the object residing in zimt:process
  // which will be passed to the act functor. This is the datum we
  // need to generate here. The second argument is a const reference
  // to the initial nD coordinate of the 'run'. This is a scalar
  // value, whereas the target is a simdized value. All we need to
  // do here is to add iota() to the 'hot' component, while the
  // other components are simply assigned, resulting in components
  // which are SIMD data, but share the same value.

  void init ( value_v & trg , const crd_t & crd ) const
  {
    trg = crd ;
    trg [ d ] += value_ele_v::iota() ;
  }

  // 'capped' variant. This is only needed if the current segment is
  // so short that no vectors can be formed at all. We fill up the
  // target value with the last valid datum, using the 'stuff'
  // member function provided for xel of simdized data.
  // If zimt::process has at least one full vector in it's current
  // 'run', the target datum (the intended input to the 'act' functor)
  // already has data which are deemed to be valid, so this (capped)
  // variant won't be called, avoiding the overhead needed for
  // 'stuffing' the unused lanes with the last 'genuine' value.
  // Why not leave the 'stuffing' to the caller (zimt::process)?
  // Because some get_t classes may, e.g., access memory to obtain
  // input, and they must be aware of the cap, to avoid 'overshooting'
  // which might result in a memory fault.

  void init ( value_v & trg ,
              const crd_t & crd ,
              std::size_t cap ) const
  {
    init ( trg , crd ) ;
    trg.stuff ( cap ) ;
  }

  // increase modifies it's argument to contain the next value.
  // The argument is, again, the 'target' datum, the intended input
  // to the 'act' functor. It's a simdized value residing in
  // zimt::process. Here, we are aware of the 'hot' axis (we received
  // it in the c'tor), so all we need to do is to increase the
  // 'hot' component by the lane count. Note how in this get_t
  // class, we don't hold any mutable state at all: the init and
  // increase member functions directly modify the target datum.
  // This is the reason for handling all arguments to these member
  // functions by reference. You may also have noticed that all
  // member functions are declared const, showing that they don't
  // change the get_t object itself.

  void increase ( value_v & trg ) const
  {
    trg [ d ] += L ;
  }

  // 'capped' variant. This is called after all vectors in the current
  // segment have been processed, so the lanes in trg beyond the cap
  // should hold valid data, and 'stuffing' them with the last datum
  // before the cap is optional. Note how we don't call 'plain'
  // increase: we want to avoid touching the values at and beyond
  // the cap, unless _stuff is true: then we fill them with the last
  // 'genuine' value.

  void increase ( value_v & trg ,
                  std::size_t cap ,
                  bool _stuff = true ) const
  {
    auto mask = ( value_ele_v::IndexesFromZero() < int(cap) ) ;
    trg [ d ] ( mask ) += T(L) ;
    if ( _stuff )
    {
      trg [ d ] ( ! mask ) = trg [ d ] [ cap - 1 ] ;
    }
  }
} ;

// class loader is a get_t which loads data from a zimt::view.
// With this get_t, 'process' can be used to implement 'coupled_f'
// We have two variants which differ only in the use of the argument
// 'stride' when calling 'bunch'. TODO: may not make a performance
// difference to call bunch with stride==1 vs. call without stride,
// in which case the unstrided variant would be superfluous.

template < typename T ,
           std::size_t N ,
           std::size_t D ,
           std::size_t L = vector_traits < T > :: vsize >
struct loader
{
  typedef zimt::xel_t < T , N > value_t ;
  typedef simdized_type < value_t , L > value_v ;
  typedef typename value_v::value_type value_ele_v ;
  typedef zimt::xel_t < long , D > crd_t ;
  typedef simdized_type < long , L > crd_v ;

  const std::size_t d ;
  const zimt::view_t < D , value_t > src ;
  const value_t * p_src ;
  const std::size_t stride ;

  // get_t's c'tor receives the zimt::view providing data and the
  // 'hot' axis. It extracts the strides from the source view.

  loader ( const zimt::view_t < D , value_t > & _src ,
           const zimt::bill_t & bill = zimt::bill_t() )
  : src ( _src ) ,
    d ( bill.axis ) ,
    stride ( _src.strides [ bill.axis ] )
  { }

  // c'tor overload for N==1. Here we also accept a view_t
  // of plain T, on top of accepting a view_t of xel_t<T,1>
  // with the general c'tor above

  template < typename = std::enable_if < N == 1 > >
  loader ( const zimt::view_t < D , T > & _src ,
           const zimt::bill_t & bill = zimt::bill_t() )
  : src ( reinterpret_cast
           < const zimt::view_t < D , value_t > & > ( _src ) ) ,
    d ( bill.axis ) ,
    stride ( _src.strides [ bill.axis ] )
  { }

  // init is used to initialize the 'target'' value to the value
  // it should hold at the beginning of the run. The discrete
  // coordinate 'crd' gives the location of the first value, and the
  // 'loader' fills it's target (the intended input for the 'act'
  // functor) with the first set of values from the source view.
  // Note how 'loader' has mutable state: the pointer p_src points
  // to data in the source view and it varies with the progess along
  // the current 'run'.

  void init ( value_v & trg , const crd_t & crd )
  {
    p_src = & ( src [ crd ] ) ;
    trg.bunch ( p_src , stride ) ;
    p_src += L * stride ;
  }

  // the 'capped' variant of 'init' only fills in the values below
  // the cap. For class loader, this is very important: there must
  // not be any memory access beyond the limit given by 'cap' to
  // avoid a memory fault.

  void init ( value_v & trg ,
              const crd_t & crd ,
              std::size_t cap )
  {
    p_src = & ( src [ crd ] ) ;
    trg.bunch ( p_src , stride , cap , true ) ;
  }

  // 'increase' fetches a full vector of data from the source view
  // and increments the pointer to the data in the view to the next
  // position.

  void increase ( value_v & trg )
  {
    trg.bunch ( p_src , stride ) ;
    p_src += L * stride ;
  }

  // The 'capped' variant only fetches a part of a full vector
  // and doesn't increase the pointer: the capped increase is
  // always the last call in a run.

  void increase ( value_v & trg ,
                  std::size_t cap ,
                  bool stuff = true )
  {
    trg.bunch ( p_src , stride , cap , stuff ) ;
  }
} ;

// unstrided_loader is a variant of class loader where the stride
// of the source view along the 'hot' axis is one.
// TODO: it may not be necessary to have this variant.
// test performance against using plain loader.

template < typename T ,
           std::size_t N ,
           std::size_t D ,
           std::size_t L = vector_traits < T > :: vsize >
struct unstrided_loader
: public loader < T , N , D , L >
{
  typedef loader < T , N , D , L > base_t ;

  using typename base_t::crd_t ;
  using typename base_t::value_t ;
  using typename base_t::value_v ;
  using base_t::src ;
  using base_t::p_src ;

  unstrided_loader ( const zimt::view_t < D , value_t > & src ,
                     const bill_t & bill = zimt::bill_t() )
  : base_t ( src , bill )
  {
    assert ( src.strides [ bill.axis ] == 1 ) ;
  }

  template < typename = std::enable_if < N == 1 > >
  unstrided_loader ( const zimt::view_t < D , T > & src ,
                     const bill_t & bill  )
  : base_t ( src , bill )
  {
    assert ( src.strides [ bill.axis ] == 1 ) ;
  }

  void init ( value_v & trg , const crd_t & crd )
  {
    p_src = & ( src [ crd ] ) ;
    trg.bunch ( p_src ) ;
    p_src += L ;
  }

  void init ( value_v & trg ,
              const crd_t & crd ,
              std::size_t cap )
  {
    p_src = & ( src [ crd ] ) ;
    trg.bunch ( p_src , 1 , cap , true ) ;
  }

  void increase ( value_v & trg )
  {
    trg.bunch ( p_src ) ;
    p_src += L ;
  }

  void increase ( value_v & trg ,
                  std::size_t cap ,
                  bool stuff = true )
  {
    trg.bunch ( p_src , 1 , cap , stuff ) ;
  }
} ;

// vloader loads vectorized data from a zimt::view_t of fundamental
// values (T). This is a good option for storing intermediate results,
// because it can use efficient SIMD load operations rather than
// having to deinterleave the data from memory holding xel of T.
// To store data so they can be loaded from vectorized storage,
// use class vstorer (see put_t.h)
// The design is so that vloader can operate within the logic used
// by zimt::process: the discrete coordinate taken by 'init' refers
// to 'where the datum would be in an array of value_t'. The array
// of T which the data are stored to has to be suitably sized;
// this is best done with the factory function zimt::get_vector_buffer.
// The array produced with this factory function has an added dimension
// zero of extent N*L, precisely of the length needed to store one
// simdized datum (of type 'value_v') along the 'hot' axis, the size
// should be so that, multiplied with the lane count, it will be
// greater or equal the 'notional' size of the unvectorized array
// - for efficiency reasons, class vstorer will only load 'full'
// simdized values from the array of simdized data. Along the other
// axes, the size must be the same as the 'notional' size of the
// equivalent unvectorized array.
// TODO: can this be 'bent' to serve as storage for tiles?

template < typename T ,
           std::size_t N ,
           std::size_t D ,
           std::size_t L = vector_traits < T > :: vsize >
struct vloader
{
  typedef zimt::xel_t < T , N > value_t ;
  typedef simdized_type < value_t , L > value_v ;

  // type of coordinate passed by the caller (zimt::process)

  typedef zimt::xel_t < long , D > crd_t ;

  // axis of the storage array which corresponds to the 'hot' axis
  // of the 'notional' array (so, d + 1 ), we'll refer to this axis
  // as the 'hot' axis of the storage array as well.

  const std::size_t d ;

  // the array serving as storage space for the vectorized data,
  // with the additional dimension zero with extent N * L

  const zimt::view_t < D + 1 , T > & src ;

  // stride along the 'hot' axis of the storage array

  const std::size_t stride ;

  // current source of the load operation

  const T * p_src ;

  // get_t's c'tor receives the zimt::view providing data and the
  // 'hot' axis. It extracts the strides from the source view.

  // vloader's c'tor receives the source array and the processing
  // axis. Here, the source array is an array of T with an added
  // dimension zero with extent N * L. The template argument 'D'
  // is the dimensionality of the 'notional' array.

  vloader ( zimt::view_t < D + 1 , T > & _src ,
            const bill_t & bill = zimt::bill_t() )
  : src ( _src ) ,
    d ( bill.axis + 1 ) ,
    stride ( _src.strides [ bill.axis + 1 ] )
  { }


  // init is used to initialize the 'target' value to the value
  // it should hold at the beginning of the run. The discrete
  // coordinate 'crd' gives the location of the first value, and the
  // 'loader' fills it's target (the intended input for the 'act'
  // functor) with the first set of values from the source view.
  // Note how 'loader' has mutable state: the pointer p_src points
  // to data in the source view and it varies with the progess along
  // the current 'run'.

  void init ( value_v & v , const crd_t & _crd , std::size_t cap = 0 )
  {
    // calculate the D+1-dimensional coordinate into 'src'.
    // This coordinate's first component is zero.

    xel_t < std::size_t , D + 1 > crd ;
    crd [ 0 ] = 0 ;

    for ( std::size_t i = 0 ; i < D ; i++ )
    {
      // the coordinate's component along the 'hot' axis is
      // divided by the lane count, the other components remain
      // the same, but all components are 'one axis further up'.

      if ( i == ( d - 1 ) )
      {
        crd [ i + 1 ] = _crd [ i ] / L ;
      }
      else
      {
        crd [ i + 1 ] = _crd [ i ] ;
      }
    }

    // with the coordinate into 'src' we can figure out the
    // initial source address

    p_src = & ( src [ crd ] ) ;
    v.load ( p_src ) ;
    p_src += stride ;
  }

  // 'increase' fetches a full vector of data from the source view
  // and increments the pointer to the data in the view to the next
  // position.

  void increase ( value_v & v , std::size_t cap = 0 ,
                  bool stuff = false )
  {
    v.load ( p_src ) ;
    p_src += stride ;
  }
} ;

// class permute is a get_t which puts together the components of
// the values it produces from separate 1D arrays, one per axis.
// This can be used to set up grids which don't have regular steps
// from one slice to the next, like for the grid_eval function
// that uses permute.

template < typename T ,
           std::size_t N ,
           std::size_t D ,
           std::size_t L = vector_traits < T > :: vsize >
struct permute
{
  typedef zimt::xel_t < T , N > value_t ;
  typedef simdized_type < value_t , L > value_v ;
  typedef typename value_v::value_type value_ele_v ;
  typedef zimt::xel_t < long , N > crd_t ;
  typedef simdized_type < long , L > crd_v ;

  const std::size_t d ;
  const std::array < zimt::view_t < 1 , T > , N > src ;
  const T * p_src ;
  const std::size_t stride ;

  // permute's c'tor receives a std::vector of 1D zimt::views providing
  // per-component data and the 'hot' axis.

  permute ( const std::array < zimt::view_t < 1 , T > , N > & _src ,
            const bill_t & bill = zimt::bill_t() )
  : src ( _src ) ,
    d ( bill.axis ) ,
    stride ( _src [ bill.axis ] . strides [ 0 ] )
  { }

  // init is used to initialize the 'target' value to the value
  // it should hold at the beginning of the run. The discrete
  // coordinate 'crd' gives the location of the first value, and the
  // 'permute' fills it's target (the intended input for the 'act'
  // functor) with the first set of values from the source view.
  // Note how 'permute' has mutable state: the pointer p_src points
  // to data in the source view and it varies with the progess along
  // the current 'run'.

  void init ( value_v & trg , const crd_t & crd )
  {
    p_src = src [ d ] . data() + crd [ d ] ;
    for ( std::size_t ch = 0 ; ch < N ; ch++ )
    {
      if ( ch == d )
      {
        trg [ ch ] . rgather ( p_src , stride ) ;
      }
      else
      {
        trg [ ch ] = src [ ch ] [ crd [ ch ] ] ;
      }
    }
    p_src += L * stride ;
  }

  // the 'capped' variant of 'init' only fills in the values below
  // the cap. For class permute, this is very important: there must
  // not be any memory access beyond the limit given by 'cap' to
  // avoid a memory fault.

  void init ( value_v & trg ,
              const crd_t & crd ,
              std::size_t cap )
  {
    p_src = src [ d ] . data() + crd [ d ] ;
    for ( std::size_t ch = 0 ; ch < N ; ch++ )
    {
      if ( ch == d )
      {
        for ( std::size_t e = 0 ; e < cap ; e++ )
          trg [ ch ] [ e ] = p_src [ e * stride ] ;
        for ( std::size_t e = cap ; e < L ; e++ )
          trg [ ch ] [ e ] = trg [ ch ] [ cap - 1 ] ;
      }
      else
      {
        trg [ ch ] = src [ ch ] [ crd [ ch ] ] ;
      }
    }
  }

  // 'increase' fetches a full vector of data from the source view
  // and increments the pointer to the data in the view to the next
  // position.

  void increase ( value_v & trg )
  {
    trg [ d ] . rgather ( p_src , stride ) ;
    p_src += L * stride ;
  }

  // The 'capped' variant only fetches a part of a full vector
  // and doesn't increase the pointer: the capped increase is
  // always the last call in a run.

  void increase ( value_v & trg ,
                  std::size_t cap ,
                  bool stuff = true )
  {
    for ( std::size_t e = 0 ; e < cap ; e++ )
      trg [ d ] [ e ] = p_src [ e * stride ] ;
    for ( std::size_t e = cap ; e < L ; e++ )
      trg [ d ] [ e ] = trg [ d ] [ cap - 1 ] ;
  }
} ;

// join_t is a get_t which loads data from N arrays of T
// (fundamentals) into simdized packets of data suitable for
// processing with the 'act' functor. This is a common scenario
// when the channels of multi-channel raster data are stored as
// separate arrays. join_t's c'tor receives a std::vector of
// N views to per-channel arrays, which should all agree in
// shape. It loads vectors of data from all arrays in turn,
// into components of the simdized datum it provides as input
// to the act functor.
// TODO: write the corresponding put_t class

template < typename T ,     // fundamental type
           std::size_t N ,  // channel count
           std::size_t D ,  // dimensions
           std::size_t L >  // lane count
struct join_t
{
  typedef zimt::xel_t < T , N > value_t ;
  typedef simdized_type < value_t , L > value_v ;
  typedef zimt::xel_t < long , D > crd_t ;

  // const std::size_t d ; // processing axis - 0 or 1 for a 2D array

  // source of data: N D-dimensional views to the fundamental type

  typedef std::array < zimt::view_t < D , T > , N > src_t ;
  src_t src ;

  zimt::xel_t < const T * , N > pickup ; // source pointers
  zimt::xel_t < long , N > stride ;      // strides of source arrays
  const std::size_t d ;

  join_t ( const src_t & _src ,
           const bill_t & bill = zimt::bill_t() )
  : src ( _src ) ,
    d ( bill.axis )
  {
    // copy out the strides of the source arrays

    for ( int ch = 0 ; ch < N ; ch++ )
    {
      stride [ ch ] = src [ ch ] . strides [ d ] ;
    }
  }

  // init is used to initialize the vectorized value to the value
  // it should hold at the beginning of the peeling run. The discrete
  // coordinate 'crd' gives the location of the first value, and this
  // function infers the start value from it. Here we set the pickup
  // pointers and perform a 'regular gather' to get the first batch
  // of values. Note how we could specialize the code to use load
  // instructions instead of rgather if the stride is 1.

  void init ( value_v & v , const crd_t & crd )
  {
    for ( int ch = 0 ; ch < N ; ch++ )
    {
      pickup [ ch ] = & ( src [ ch ] [ crd ] ) ;
      v [ ch ] . rgather ( pickup [ ch ] , stride [ ch ] ) ;
    }
  }

  // 'capped' variant. This is only needed if the current segment is
  // so short that no vectors can be formed at all. We fill up the
  // target value with the last valid datum.

  void init ( value_v & trg ,
              const crd_t & crd ,
              std::size_t cap )
  {
    for ( int ch = 0 ; ch < N ; ch++ )
    {
      pickup [ ch ] = & ( src [ ch ] [ crd ] ) ;
      for ( std::size_t e = 0 ; e < cap ; e++ )
        trg [ ch ] [ e ] = pickup [ ch ] [ e * stride [ ch ] ] ;
      trg.stuff ( cap ) ;
    }
  }


  // increase modifies it's argument to contain the next value, or
  // next vectorized value, respectively - first we increase the
  // pickup pointers, then we get the data from that location.

  void increase ( value_v & trg )
  {
    for ( int ch = 0 ; ch < N ; ch++ )
    {
      pickup [ ch ] += L * stride [ ch ] ;
      trg [ ch ] . rgather ( pickup [ ch ] , stride [ ch ] ) ;
    }
  }

  // 'capped' variant. This is called after all vectors in the current
  // segment have been processed, so the lanes in trg beyond the cap
  // should hold valid data, and 'stuffing' them with the last datum
  // before the cap is optional.

  void increase ( value_v & trg ,
                  std::size_t cap ,
                  bool _stuff = true )
  {
    for ( int ch = 0 ; ch < N ; ch++ )
    {
      pickup [ ch ] += L * stride [ ch ] ;
      for ( std::size_t e = 0 ; e < cap ; e++ )
        trg [ ch ] [ e ] = pickup [ ch ] [ e * stride [ ch ] ] ;
      if ( _stuff )
        trg.stuff ( cap ) ;
    }
  }
} ;

// linspace_t produces input to the 'act' functor which is
// like a NumPy ngrid: The data start out at a given initial
// value and increase (or decrease) by a fixed amount along
// each axis. It's also possible to pass discrete coordinates
// on to the 'act' functor and do the range transformation
// there, but the 'act' functor should not have a notion of
// which axis is 'hot', but rather take the entire simdized
// datum it receives as input 'at face value', so doing the
// transformation there would involve all channels, whereas
// linspace_t is aware of the 'hot' axis and can produce the
// simdized data more efficiently, by only modifying the
// component of the simdized datum pertaining to the hot axis.

template < typename T ,     // fundamental type
           std::size_t N ,  // channel count
           std::size_t D ,  // dimensions
           std::size_t L >  // lane count
struct linspace_t
{
  typedef zimt::xel_t < T , N > value_t ;
  typedef simdized_type < value_t , L > value_v ;
  typedef typename value_v::value_type value_ele_v ;
  // I'd like to use long, but I don't get it vectorized everywhere
  typedef zimt::xel_t < int , D > crd_t ;
  typedef simdized_type < crd_t , L > crd_v ;
  typedef typename crd_v::value_type crd_ele_v ;

  const std::size_t d ;
  const value_t start ;
  const value_t step ;

  // linspace_t's c'tor receives start, step and axis. Note how
  // start and step are N-dimensional; each component gives the
  // intended value for the corresponding axis.

  linspace_t ( const value_t & _start ,
               const value_t & _step ,
               const bill_t & bill = zimt::bill_t() )
  : start ( _start ) ,
    step ( _step ) ,
    d ( bill.axis )
  { }

  // init is used to initialize the vectorized value to the value
  // it should hold at the beginning of the peeling run. The discrete
  // coordinate 'crd' gives the location of the first value, and this
  // function infers the start value from it. The scalar value will
  // not be used until peeling is done, so it isn't initialized here.

  void init ( value_v & trg , const crd_t & crd )
  {
    trg = step * crd + start ;
    trg [ d ] += value_ele_v::iota() * step [ d ] ;
  }

  // 'capped' variant. This is only needed if the current segment is
  // so short that no vectors can be formed at all. We fill up the
  // target value with the last valid datum.

  void init ( value_v & trg ,
              const crd_t & crd ,
              const std::size_t & cap )
  {
    init ( trg , crd ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;

    // to 'go over the lanes' instead:

    // if ( cap == L )
    // {
    //   init ( trg , crd ) ;
    // }
    // else
    // {
    //   trg = step * crd + start ;
    //   for ( std::size_t e = 0 ; e < cap ; e++ )
    //     trg [ d ] [ e ] += T ( e ) * step [ d ] ;
    //   trg.stuff ( cap ) ;
    // }
  }

  // increase modifies it's argument to contain the next value

  void increase ( value_v & trg )
  {
    trg [ d ] += ( step [ d ] * L ) ;
  }

  // 'capped' variant. This is called after all vectors in the current
  // segment have been processed. If we'd 'go over the lanes', we might
  // only set entries below cap and rely on the higher lanes holding
  // viable values from processing the last full vector, but it's
  // probably more efficient to fill the target with a vector operation
  // and stuff unconditionally.

  void increase ( value_v & trg ,
                  const std::size_t & cap ,
                  const bool & _stuff = true )
  {
    increase ( trg ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;

    // to go 'over the lanes' instead:

    // if ( cap == L )
    // {
    //   increase ( trg ) ;
    // }
    // else
    // {
    //   for ( std::size_t e = 0 ; e < cap ; e++ )
    //     trg [ d ] [ e ] += ( step [ d ] * L ) ;
    //   if ( _stuff )
    //     trg.stuff ( cap ) ;
    // }
  }
} ;

template < typename T ,     // fundamental type
           std::size_t N ,  // channel count
           std::size_t D ,  // dimensions
           std::size_t L >  // lane count
struct gridspace_t
{
  typedef zimt::xel_t < T , N > value_t ;
  typedef simdized_type < value_t , L > value_v ;
  typedef typename value_v::value_type value_ele_v ;
  typedef zimt::xel_t < long , D > crd_t ;
  typedef simdized_type < crd_t , L > crd_v ;
  typedef typename crd_v::value_type crd_ele_v ;

  const std::size_t d ;
  const value_t start ;
  const std::array < value_t , D > step ;

  // gridspace_t's c'tor receives start, step and axis. Note how
  // start is a single N-dimensional value, whereas step has one
  // value for the increase for a step in each of the canonical
  // directions. This way, it can produce grids which are not
  // perpendicular to the canonical directions.

  gridspace_t ( const value_t & _start ,
                const std::array < value_t , D > & _step ,
                const bill_t & bill = zimt::bill_t() )
  : start ( _start ) ,
    step ( _step ) ,
    d ( bill.axis )
  { }

  // init is used to initialize the vectorized value to the value
  // it should hold at the beginning of the peeling run. The discrete
  // coordinate 'crd' gives the location of the first value, and this
  // function infers the start value from it.

  void init ( value_v & trg , const crd_t & crd )
  {
    value_t initial = start ;

    // for each axis a, extract the increase vector and multiply it
    // with the coordinate's component in direction a. Add these
    // per-axis displacemets to 'initial', yielding the first lane
    // of the target vector.

    for ( int a = 0 ; a < D ; a++ )
      initial += crd [ a ] * step [ a ] ;

    // use that to initialize the entire target vector

    trg = initial ;

    // and add increments along the hot axis

    trg += value_ele_v::iota() * step [ d ] ;
  }

  // 'capped' variant. This is only needed if the current segment is
  // so short that no vectors can be formed at all. We fill up the
  // target value with the last valid datum. Since the values in trg
  // are generated by simple arithmetic, there's no harm in first
  // filling the vector 'as usual' and then overwriting the ones
  // beyond the cap with the last valid entry. Alternatively, we
  // could proceed as in the previous class, going 'over the lanes'

  void init ( value_v & trg ,
              const crd_t & crd ,
              const std::size_t & cap )
  {
    init ( trg , crd ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }

  // increase modifies it's argument to contain the next value. Note how
  // step[d] is a vector with D components, rather than the single delta
  // value it was in the previous class. The increase will affect all
  // of the target's components.

  void increase ( value_v & trg )
  {
    trg += ( step [ d ] * L ) ;
  }

  // 'capped' variant. This is called after all vectors in the current
  // segment have been processed. Again we generate the target 'as usual'
  // and then call stuff unconditionally, rather than 'going over the
  // lanes'.

  void increase ( value_v & trg ,
                  const std::size_t & cap ,
                  const bool & _stuff = true )
  {
    increase ( trg ) ;
    if ( cap < L )
      trg.stuff ( cap ) ;
  }
} ;

// no_input is a get_t which doesn't 'touch' md_in in
// zimt::process. The 'act' functor can refer to it, but it's
// content is not defined. So what's the point? The act functor
// might be made to do some task at every invocation which does
// not depend on a location, i.e. produce or process a random
// value.

template < typename T ,     // fundamental type
           std::size_t N ,  // channel count
           std::size_t D ,  // dimensions
           std::size_t L >  // lane count
struct no_input
{
  template < typename ... A >
  void init ( A ... args )
  { }

  template < typename ... A >
  void increase ( A ... args )
  { }
} ;

// to avoid having to deal with the concrete type of a get_t, we can
// use type erasure, a technique I call 'grokking'. Basically, it
// captures an object's set of member functions as std::functions
// and provides a new object which delegates to these std::functions.
// In result the new object is decoupled from the 'grokkee' type, but
// provides it's functionality just the same. The internal workings
// do still use the 'grokked' object, but it's type is hidden from
// view - it's been 'erased'.
// For a get_t, we need an object providing two init and two increase
// overloads - both for the 'normal' and the capped variants of init
// and increase. Here's the class definition - it's quite a mouthful,
// but further down there's a factory function to perform the 'grok'
// which uses ATD, making the process simple.

template < typename T ,     // fundamental type
           std::size_t N ,  // channel count
           std::size_t D ,  // dimensions
           std::size_t L >  // lane count
class grok_get_t
: private grok_t
{
  // we need some of the grokkee's types

  typedef zimt::xel_t < T , N > value_t ;
  typedef simdized_type < value_t , L > value_v ;
  typedef zimt::xel_t < long , D > crd_t ;

  // grok_get_t holds four std::functions:

  std::function < void ( void * & ,
                         value_v & ,
                         const crd_t & ) > vinit ;

  std::function < void ( void * & ,
                         value_v & ) > vincrease ;

  std::function < void ( void * & ,
                         value_v & ,
                         const crd_t & ,
                         const std::size_t & ) > cinit ;

  std::function < void ( void * & ,
                         value_v & ,
                         const std::size_t & ,
                         const bool & ) > cincrease ;

public:

  template < typename grokkee_t >
  grok_get_t ( const grokkee_t & grokkee )
  : grok_t ( grokkee )
  {
    // the std::functions are initialized with wrappers taking
    // p_context and a set of arguments which are passed on to
    // the grokkee's member functions.

    vinit = [] ( void * & p_ctx ,
                 value_v & v ,
                 const crd_t & crd )
          {
            auto p_gk = static_cast<grokkee_t*> ( p_ctx ) ;
            p_gk->init ( v , crd ) ;
          } ;

    cinit = [] ( void * & p_ctx ,
                 value_v & v ,
                 const crd_t & crd ,
                 const std::size_t & cap )
          {
            auto p_gk = static_cast<grokkee_t*> ( p_ctx ) ;
            p_gk->init ( v , crd , cap ) ;
          } ;

    vincrease = [] ( void * & p_ctx ,
                     value_v & v )
          {
            auto p_gk = static_cast<grokkee_t*> ( p_ctx ) ;
            p_gk->increase ( v ) ;
          } ;

    cincrease = [] ( void * & p_ctx ,
                     value_v & v ,
                     const std::size_t & cap ,
                     const bool & stuff )
          {
            auto p_gk = static_cast<grokkee_t*> ( p_ctx ) ;
            p_gk->increase ( v , cap , stuff ) ;
          } ;
  }

  // grok_get_t itself offers the typical member functions of
  // a get_t object, which in turn delegate to the stored
  // std::functions.

  void init ( value_v & trg , const crd_t & crd )
  {
    vinit ( p_context , trg , crd ) ;
  }

  void init ( value_v & trg ,
              const crd_t & crd ,
              const std::size_t & cap )
  {
    cinit ( p_context , trg , crd , cap ) ;
  }

  void increase ( value_v & trg )
  {
    vincrease ( p_context , trg ) ;
  }

  void increase ( value_v & trg ,
                  const std::size_t & cap ,
                  const bool & _stuff = true )
  {
    cincrease ( p_context , trg , cap , _stuff ) ;
  }
} ;

// grok_get is a factory function to 'grok' a get_t object.
// using ATD, the invocation to 'grok' some get_t x is simply
// auto gk = grok_get ( x ) ;
// gk can then be used wherever a get_t is required.

template < typename T ,     // fundamental type
           std::size_t N ,  // channel count
           std::size_t D ,  // dimensions
           std::size_t L ,  // lane count
           template < typename ,
                      std::size_t ,
                      std::size_t ,
                      std::size_t >
             class G >
grok_get_t < T , N , D , L > grok_get
  ( G < T , N , D , L > grokkee )
{
  return grok_get_t < T , N , D , L > ( grokkee ) ;
}

// 'suffixed' is a composite of a get_t and an act-like functor
// post-processing the output of the get_t. The resulting object
// behaves as a get_t.

template < typename T ,     // fundamental type
           std::size_t N ,  // channel count
           std::size_t D ,  // dimensions
           std::size_t L >  // lane count
struct suffixed_t
{
  typedef zimt::xel_t < long , D > crd_t ;
  typedef zimt::xel_t < T , N > value_t ;
  typedef simdized_type < value_t , L > value_v ;

  grok_get_t < T , N , D , L > head ;
  grok_type < value_t , value_t , L > tail ;

  suffixed_t ( const grok_get_t < T , N , D , L > & _head ,
               const grok_type < value_t , value_t , L > & _tail )
  : head ( _head ) ,
    tail ( _tail )
  { }

  void init ( value_v & trg , const crd_t & crd )
  {
    head.init ( trg , crd ) ;
    tail.eval ( trg , trg ) ;
  }

  void init ( value_v & trg ,
              const crd_t & crd ,
              const std::size_t & cap )
  {
    head.init ( trg , crd , cap ) ;
    tail.eval ( trg , trg ) ;
  }

  void increase ( value_v & trg )
  {
    head.increase ( trg ) ;
    tail.eval ( trg , trg ) ;
  }

  void increase ( value_v & trg ,
                  const std::size_t & cap ,
                  const bool & _stuff = true )
  {
    head.increase ( trg , cap , _stuff ) ;
    tail.eval ( trg , trg ) ;
  }
} ;

// overload which takes grokkee type with different type signature.
// to use this overload, you'll have to pass template arguments
// T, N, D and L.

template < typename T ,     // fundamental type
           std::size_t N ,  // channel count
           std::size_t D ,  // dimensions
           std::size_t L ,  // lane count
           class G >
grok_get_t < T , N , D , L > grok_get ( G grokkee )
{
  return grok_get_t < T , N , D , L > ( grokkee ) ;
}

// fusion_t is used to form a synopsis of the values produced by
// several get_t objects. The get_t objects are passed to fusion_t's
// c'tor as a std::vector of grok_get_t, so the get_t objects which are
// put to use can be of different types, but their argument types have
// to agree, so that they can be held by the same type of grok_get_t.
// fusion_t itself acts as a get_t. All of it's member functions
// invoke the grok_get_t's corresponding member functions, and the
// results are stored in a std::vector of the grok_get_t's result type
// (a simdized datum, here called 'partial_v'). Once all the grok_get_t
// have been invoked and src_v - the std::vector of partial_v - is filled,
// the 'synopsis' callback is called. It receives a const& to src_v,
// and a reference to the final result of the fusion_t's operation -
// the datum which is passed on to the 'act' functor in the zimt::process
// invocation. The synopsis callback handles the reduction from the set
// of results of the 'inner' grok_get_t objects into the 'outer' result.
// An example: if the grok_get_t objects in 'get_v' produce coordinates,
// the synopsis function might have a set of interpolators to produce
// pixel values from coordinates. It might feed each interpolator with
// a corresponding coordinate and form a weighted sum of the results,
// which it writes to the final result - the synopsis. This is received
// by the act functor and processing continues as usual. An alternative
// modus operandi for the synopsis function would be a 'winner takes all'
// scheme, where the synopsis function picks the 'best' one from the
// set of values produced by the grok_get_t and only evaluates the
// single corresponding interpolator to produce a pixel output.
// Note how - in this example of fusing interpolated values - the value_v
// held in zimt::process is not of the type of the coordinates (that's
// partial_t/partial_v) but instead value_t/value_v which is the result
// of the synopsis function, so when this is part of a pixel pipeline,
// the act functor in zimt::process needs to be set up to accept a
// pixel value rather than a coordinate value - it may be sufficient to
// use a pass_through, if the pixel is already 'ready' at this stage.
// The SYN-type object is the synopsis-formin agent. Initially I coded
// this functionality as a callback function, but I think its' better
// to set this up as part of this object. The SYN-type object must
// provide operator() overloads for each type of partial_v which it
// can handle, with this signature:
// void operator() ( const std::vector < partial_v > & ,
//                   value_v & trg ,
//                   const std::size_t & )

template < typename T ,     // fundamental type of result (value_t)
           std::size_t N ,  // channel count of result
           std::size_t D ,  // dimensions
           std::size_t L ,  // lane count
           typename U ,     // fundamental type of partial_t
           std::size_t M ,  // channel count of partial_t
           typename SYN >   // class of synopsis-forming object
struct fusion_t
{
  typedef zimt::xel_t < T , N > value_t ;
  typedef simdized_type < value_t , L > value_v ;
  typedef zimt::xel_t < long , D > crd_t ;
  typedef zimt::xel_t < U , M > partial_t ;
  typedef simdized_type < partial_t , L > partial_v ;
  typedef grok_get_t < U , M , D , L > gg_t ;
  SYN synopsis ;

  std::vector < gg_t > get_v ;
  std::vector < partial_v > src_v ;
  const std::size_t sz ;

  fusion_t ( const std::vector < gg_t > & _get_v ,
             const SYN & _synopsis )
  : src_v ( _get_v.size() ) ,
    get_v ( _get_v ) ,
    sz ( _get_v.size() ) ,
    synopsis ( _synopsis )
  { }

  void init ( value_v & trg , const crd_t & crd )
  {
    for ( std::size_t i = 0 ; i < sz ; i++ )
    {
      get_v[i].init ( src_v[i] , crd ) ;
    }
    synopsis ( src_v , trg , L ) ;
  }

  void init ( value_v & trg ,
              const crd_t & crd ,
              std::size_t cap )
  {
    for ( std::size_t i = 0 ; i < sz ; i++ )
    {
      get_v[i].init ( src_v[i] , crd , cap ) ;
    }
    synopsis ( src_v , trg , cap ) ;
  }

  void increase ( value_v & trg )
  {
    for ( std::size_t i = 0 ; i < sz ; i++ )
    {
      get_v[i].increase ( src_v[i] ) ;
    }
    synopsis ( src_v , trg , L ) ;
  }

  void increase ( value_v & trg ,
                  std::size_t cap ,
                  bool _stuff = true )
  {
    for ( std::size_t i = 0 ; i < sz ; i++ )
    {
      get_v[i].increase ( src_v[i] , cap , _stuff ) ;
    }
    synopsis ( src_v , trg , cap ) ;
  }
} ;

// class zip_t 'harnesses' two get_t in sync to produce their
// respective results. These results are passed to the synopsis
// functor, which in turn produces a result which constitutes
// the output of the zip_t object. Thy synopsis functor only
// needs to provide a single (vector) function - a specific
// type is not enforced, and a good choice is a lambda. An
// alternative would be to use a std::function which would
// allow for more restrictive parameterization.
// The first four template arguments are 'standard fare',
// followed by the synopsis object's type (use decltype to
// specify a lambda's type). The remaining template arguments
// fix fundamental type and channel count of the second get_t
// and result, in case they differ.
// Note that the SYN object should provide a call signature
// with a trailing 'cap' argument with a default value, like:
//   auto syn = [] ( const x_v & v1 , const y_v & v2 ,
//                   z_v & v3 , std::size_t cap = LANES )

template < typename T1 ,    // fundamental type of first get_t
           std::size_t N1 , // channel count of of first get_t
           std::size_t D ,  // dimensions
           std::size_t L ,  // lane count
           typename SYN ,   // class of synopsis-forming object
           typename T2 = T1 , // fundamental type of second get_t
           std::size_t N2 = N1 , // channel count of second get_t
           typename T3 = T1 , // fundamental type of result
           std::size_t N3 = N1 > // channel count of result
struct zip_t
{
  typedef zimt::xel_t < T1 , N1 > value1_t ;
  typedef simdized_type < value1_t , L > value1_v ;
  typedef zimt::xel_t < T2 , N2 > value2_t ;
  typedef simdized_type < value2_t , L > value2_v ;
  typedef zimt::xel_t < T3 , N3 > value3_t ;
  typedef simdized_type < value3_t , L > value3_v ;

  typedef zimt::xel_t < long , D > crd_t ;

  typedef grok_get_t < T1 , N1 , D , L > gg1_t ;
  typedef grok_get_t < T2 , N2 , D , L > gg2_t ;
  
  gg1_t get1 ;
  gg2_t get2 ;
  SYN synopsis ;

  zip_t ( const gg1_t & _get1 ,
          const gg2_t & _get2 ,
          const SYN & _synopsis )
  : get1 ( _get1 ) ,
    get2 ( _get2 ) ,
    synopsis ( _synopsis )
  { }

  void init ( value3_v & trg , const crd_t & crd )
  {
    value1_v v1 ;
    value2_v v2 ;
    get1.init ( v1 , crd ) ;
    get2.init ( v2 , crd ) ;
    synopsis ( v1 , v2 , trg ) ;
  }

  void init ( value3_v & trg ,
              const crd_t & crd ,
              std::size_t cap )
  {
    value1_v v1 ;
    value2_v v2 ;
    get1.init ( v1 , crd , cap ) ;
    get2.init ( v2 , crd , cap ) ;
    synopsis ( v1 , v2 , trg , cap ) ;
  }

  void increase ( value3_v & trg )
  {
    value1_v v1 ;
    value2_v v2 ;
    get1.increase ( v1 ) ;
    get2.increase ( v2 ) ;
    synopsis ( v1 , v2 , trg ) ;
  }

  void increase ( value3_v & trg ,
                  std::size_t cap ,
                  bool _stuff = true )
  {
    value1_v v1 ;
    value2_v v2 ;
    get1.increase ( v1 , cap , _stuff ) ;
    get2.increase ( v2 , cap , _stuff ) ;
    synopsis ( v1 , v2 , trg , cap ) ;
  }
} ;

END_ZIMT_SIMD_NAMESPACE
HWY_AFTER_NAMESPACE() ;

#endif // sentinel
