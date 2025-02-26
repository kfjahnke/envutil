/************************************************************************/
/*                                                                      */
/*   utility to convert and extract images from 360 degree environments */
/*                                                                      */
/*            Copyright 2024 by Kay F. Jahnke                           */
/*                                                                      */
/*    The git repository for this software is at                        */
/*                                                                      */
/*    https://github.com/kfjahnke/envutil                               */
/*                                                                      */
/*    Please direct questions, bug reports, and contributions to        */
/*                                                                      */
/*    kfjahnke+envutil@gmail.com                                        */
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

// this header provides a parser for PTO files. The PTO syntax is quite
// simple, so the parser only uses two regular expressions: one for entire
// lines in the PTO file, and one for the fields inside each line. The
// parser only processes lines which start with a letter - other lines
// are comments, starting with '#', or empty, or don't start with a letter.
// the fields always start with one or several letters and are followed
// immediately by a value, which may be in quotes for string values.

#include <map>
#include <vector>
#include <fstream>
#include <regex>

struct pto_line_type
{
  std::string original ;
  std::string head ; // will contain the 'header', like "i" for image lines
  std::map < std::string , std::string > field_map ; // map of items in the line
} ;

struct pto_parser_type
{
  const std::regex pto_line_regex ; // matches a whole pto line
  const std::regex pto_item_regex ; // matches a single item in a line

  // holds groups of lines starting with the same letter

  std::map < std::string , std::vector < pto_line_type > > line_group ;

  // the c'tor initializes the REs

  pto_parser_type()
  : pto_line_regex ( "([a-zA-Z])\\s(.+)[\n\r]*" ) ,
    pto_item_regex ( "([A-Za-z]+)((\"[^\"]+\")|(\\S*))" )
  { }

  // parse_pto_line receives a line from a PTO file and scans it's content.
  // the result is saved as a pto_line_type object. This is merely the
  // textual representation of the PTO file, exttracting the values is
  // done in a second step.
  // The data are stored in groups: all lines headed by the same letter
  // are stored in a std::vector of pto_line, which itself is stored in
  // a std::map indexed with that heading letter. This scheme preserves
  // the order in which the lines are found, which is relevant because
  // there are sometimes references to previous lines.

  bool parse_pto_line ( const std::string & s )
  {
    // split the line into head and remainder

    std::smatch parts ;
    auto success = std::regex_match ( s , parts , pto_line_regex ) ;
    if ( ! success )
    {
      std::cout << "parser ignores: " << s << std::endl ;
      return true ;
    }

    pto_line_type pto_line ;

    std::ssub_match part = parts[1] ;
    pto_line.head = part.str() ;
    pto_line.original = s ;
    std::string tail = parts[2].str() ;

    // now split the remainder into individual items

    auto start = std::sregex_iterator ( tail.begin() ,
                                        tail.end() ,
                                        pto_item_regex ) ;
    auto end = std::sregex_iterator() ;

    for ( auto i = start ; i != end ; ++i )
    {
      // iterate over the items and separate name and value
  
      // note: since the regex_iterator uses regex_search, all items
      // we encounter now are RE matches - if the input is malformed,
      // it's skipped over. This may be an issue.

      auto item = i->str() ;
      std::regex_match ( item , parts , pto_item_regex ) ;
      auto field_name = parts[1].str() ;
      auto field_value = parts[2].str() ;

      // store the field, value pair in the field map

      if ( field_value[0] == '=' )
      {
        auto number = field_value.substr ( 1 , field_value.size() - 1 ) ;
        int refers_to = std::stoi ( number ) ;
        if ( field_name != "j" )
        {
          auto source_line = line_group [ "i" ] [ refers_to ] ;
          field_value = source_line.field_map [ field_name ] ;
        }
        pto_line.field_map [ field_name ] = field_value ;
      }
      else
        pto_line.field_map [ field_name ] = field_value ;
    }

    auto it = line_group.find ( pto_line.head ) ;
    if ( it == line_group.end() )
    {
      line_group [ pto_line.head ] = std::vector < pto_line_type >() ;
      it = line_group.find ( pto_line.head ) ;
    }
    it->second.push_back ( pto_line ) ;
    return true ;
  }

  // TODO: handle escaped quotation marks

  bool read_pto_file ( const std::string & filename )
  {
    std::ifstream str ( filename ) ;
    if ( ! str )
    {
      std::cerr << "could not open pto file " << filename << std::endl ;
      return false ;
    }

    // ought to suffice

    char buffer [ 2048 ] ;

    while ( str.getline ( buffer , 2048 ) )
    {
      bool success = parse_pto_line ( buffer ) ;
      if ( ! success )
        return false ;
    }

    return true ;
  }

  void walk()
  {
    for ( const auto & group_it : line_group )
    {
      std::cout << "line group: " << group_it.first << std::endl ;
      const auto & line_list ( group_it.second ) ;
      int lineno = 0 ;
      
      for ( const auto & line_it : line_list )
      {
        std::cout << lineno << ": " << line_it.original << std::endl ;
        ++ lineno ;
        for ( const auto & field_it: line_it.field_map )
        {
          auto content = field_it.second ;
          if ( content[0] == '=' )
          {
            auto number = content.substr ( 1 , content.size() - 1 ) ;
            int refers_to = std::stoi ( number ) ;
            if ( field_it.first == "j" )
            {
              std::cout << "    " << field_it.first
                       << " " << refers_to
                        << std::endl ;
            }
            else
            {
              auto source_line = line_list [ refers_to ] ;
              std::cout << "    " << field_it.first
                        << " " << source_line.field_map [ field_it.first ]
                        << std::endl ;
            }
          }
          else
          {
            std::cout << "    " << field_it.first
                      << " " << field_it.second << std::endl ;
          }
        }
      }
    }
  }
} ;

