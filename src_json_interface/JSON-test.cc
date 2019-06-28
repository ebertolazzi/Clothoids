//
//  main.cpp
//  genconjson
//
//  Created by Nicola Dal Bianco on 08/11/17.
//  Copyright © 2017 Nicola Dal Bianco. All rights reserved.
//

#include "GenericContainerJsonHandler.hh"
#include <iostream>
#include <string>

using namespace GC;
using namespace std;

int main( int argc, const char * argv[] )
{

  GenericContainer gc;

  int test = 1; // 1 or 2
  switch ( test ) {
  case 1: {
    gc["a"] = 2;
    gc["b"] = vector<int>( {1, 2, 3, 4, 5} );
    gc["bb"].set_vec_string();
    gc["bb"] = vector<real_type>( {4.5, 0.000000000000000000009086954454768098980679} );
    gc["c"] = false;
    GenericContainer gc2;
    gc2["asd"].set_vector();
    gc2["asd"].get_vector().push_back( GenericContainer() );
    gc2["asd"].get_vector().push_back( GenericContainer() );
    gc2["asd"][0].set_vec_int( {0, 1, 2, 4} );
    gc2["asd"][1].set_vec_string( {"1+2i", "3.4", "34i"} );
    gc["d"] = gc2;
    gc2["eee"] = vector<bool>( {false, true} );
    mat_real_type mat( 2, 2 );
    mat( 0, 0 ) = 0;
    mat( 1, 0 ) = 1;
    mat( 0, 1 ) = 2.2;
    mat( 1, 1 ) = 3;
    gc["e"].set_mat_real( mat );
    gc["f"];
    GenericContainer ff;
    ff.set_vector();
    ff.push_int( 3 );
    ff.push_bool( true );
    gc["g"] = ff;
    break;
  }
  case 2: {
    GC::vector_type & v = gc.set_vector();
    v.resize( 10 );
    v[0] = 1;
    v[1].set_vec_real();
    v[2].set_map();
    v[3].set_vec_string();
    v[4] = 1.3;
    v[5] = "pippo";
    v[6].set_map();
    v[7].set_vector();
    v[8] = true;
    GC::vec_real_type & vv = v[1].get_vec_real();
    vv.resize( 10 );
    vv[2] = 123;
    GC::map_type & mm = v[2].get_map();
    mm["pippo"]    = 13;
    mm["pluto"]    = 1;
    mm["paperino"] = 3;
    GenericContainer & gmm = v[2]; // access element 2 as GenericContainer
    gmm["aaa"]     = "stringa1";   // is the same as mm["aaa"] = "stringa"
    gmm["bbb"]     = "stringa2";   // is the same as mm["aaa"] = "stringa"
    GC::vec_string_type & vs = v[3].get_vec_string();
    vs.push_back( "string1" );
    vs.push_back( "string2" );
    vs.push_back( "string3" );
    vs.push_back( "string4" );
    GC::map_type & m = v[6].get_map();
    m["aaa"]    = 123;
    m["bbb"]    = 3.4;
    m["vector"].set_vec_int();
    GC::vec_int_type & vi = m["vector"].get_vec_int();
    vi.push_back( 12 );
    vi.push_back( 10 );
    vi.push_back( 1 );

    GC::vector_type & vg = v[7].get_vector();
    vg.resize( 4 );
    vg[0] = 123;
    vg[1] = 3.14;
    vg[2] = "nonna papera";
    break;
  }
  default:
    break;
  }

  cout << "Original container is:" << endl;
  gc.dump( cout );
  cout << endl << endl;

  GenericContainer gc_options;
  //    gc_options[GC_JSON_PRETTY] = true;
  //    gc_options[GC_JSON_INDENT_CHAR] = "\t";
  //    gc_options[GC_JSON_INDENT_NUM] = 2;
  //    gc_options[GC_JSON_MAT_ORDER] = row_major;
  string out_str = genericContainerToJsonString( gc, gc_options );
  cout << "Json string is:\n" << out_str << std::endl;
  cout << endl << endl;

  GenericContainer gc_back;
  jsonStringToGenericContainer( out_str, gc_back, gc_options );

  cout << "Re-converted container is (note complex string are converted to complex numbers):" << endl;
  gc_back.print( cout );
  cout << endl << endl;

  // test a generic json string
  string a_json = "[  {    \"_id\": \"5a0ae30d8419c393f642fec3\",    \"index\": 0,    \"guid\": \"d0268ff7-cde2-48e4-9c75-5740c53c3b23\",    \"isActive\": true,    \"balance\": \"$1,949.87\",    \"picture\": \"http://placehold.it/32x32\",    \"age\": 26,    \"eyeColor\": \"green\",    \"name\": \"Sargent Frazier\",    \"gender\": \"male\",    \"company\": \"MANGLO\",    \"email\": \"sargentfrazier@manglo.com\",    \"phone\": \"+1 (894) 561-3896\",    \"address\": \"769 Prospect Avenue, Snyderville, North Dakota, 3198\",    \"about\": \"Qui ad eu officia commodo aute dolor proident ea esse fugiat deserunt sint anim incididunt. Cupidatat ipsum tempor ipsum cillum laborum eu culpa eiusmod laborum quis irure in proident. Amet ipsum aute nulla et. Ipsum sit velit cillum in consequat. Est exercitation deserunt ad Lorem laborum occaecat mollit cupidatat fugiat quis sunt elit voluptate.\",    \"registered\": \"2016-01-07T03:21:06 -01:00\",    \"latitude\": -38.243575,    \"longitude\": -0.862908,    \"tags\": [      \"sunt\",      \"voluptate\",      \"labore\",      \"sunt\",      \"elit\",      \"cupidatat\",      \"nulla\"    ],    \"friends\": [      {        \"id\": 0,        \"name\": \"April Dunlap\"      },      {        \"id\": 1,        \"name\": \"Deirdre Mayer\"      },      {        \"id\": 2,        \"name\": \"Bernadette Durham\"      }    ],    \"greeting\": \"Hello, Sargent Frazier! You have 7 unread messages.\",    \"favoriteFruit\": \"strawberry\"  },  {    \"_id\": \"5a0ae30da2e6848d6f374bba\",    \"index\": 1,    \"guid\": \"578ab333-ea26-40c2-a130-04d47aa9000d\",    \"isActive\": true,    \"balance\": \"$2,690.13\",    \"picture\": \"http://placehold.it/32x32\",    \"age\": 26,    \"eyeColor\": \"blue\",    \"name\": \"Vickie Roy\",    \"gender\": \"female\",    \"company\": \"INSOURCE\",    \"email\": \"vickieroy@insource.com\",    \"phone\": \"+1 (996) 421-2072\",    \"address\": \"922 Dahill Road, Byrnedale, Marshall Islands, 1410\",    \"about\": \"Adipisicing amet veniam dolore in exercitation consequat consectetur cupidatat enim non. Sint minim velit deserunt in aliqua excepteur nostrud ea cupidatat fugiat excepteur mollit consequat. Tempor do minim qui qui labore cupidatat.\",    \"registered\": \"2016-12-20T05:28:05 -01:00\",    \"latitude\": -13.097375,    \"longitude\": -144.934101,    \"tags\": [      \"quis\",      \"duis\",      \"do\",      \"magna\",      \"qui\",      \"deserunt\",      \"nostrud\"    ],    \"friends\": [      {        \"id\": 0,        \"name\": \"Riddle Burt\"      },      {        \"id\": 1,        \"name\": \"Hughes Howell\"      },      {        \"id\": 2,        \"name\": \"Mcbride Tyson\"      }    ],    \"greeting\": \"Hello, Vickie Roy! You have 3 unread messages.\",    \"favoriteFruit\": \"banana\"  },  {    \"_id\": \"5a0ae30df49a05d7a259f116\",    \"index\": 2,    \"guid\": \"84f91917-7331-441d-aa3c-da08994bc1a9\",    \"isActive\": true,    \"balance\": \"$3,503.13\",    \"picture\": \"http://placehold.it/32x32\",    \"age\": 25,    \"eyeColor\": \"brown\",    \"name\": \"Delacruz Schultz\",    \"gender\": \"male\",    \"company\": \"ISOLOGIA\",    \"email\": \"delacruzschultz@isologia.com\",    \"phone\": \"+1 (949) 553-3319\",    \"address\": \"240 Homecrest Avenue, Groveville, Alaska, 3106\",    \"about\": \"Eu eiusmod nulla incididunt eu mollit et id elit nisi velit Lorem. Id Lorem in fugiat proident dolor. Exercitation elit Lorem ipsum pariatur duis do consectetur ex laboris incididunt. Exercitation nisi sit ea ea ullamco. Excepteur ipsum duis et id excepteur.\",    \"registered\": \"2016-02-14T03:02:53 -01:00\",    \"latitude\": -7.842765,    \"longitude\": -166.808559,    \"tags\": [      \"velit\",      \"esse\",      \"ad\",      \"laborum\",      \"incididunt\",      \"ex\",      \"fugiat\"    ],    \"friends\": [      {        \"id\": 0,        \"name\": \"Romero Hodge\"      },      {        \"id\": 1,        \"name\": \"Dona Reilly\"      },      {        \"id\": 2,        \"name\": \"Earlene Richards\"      }    ],    \"greeting\": \"Hello, Delacruz Schultz! You have 9 unread messages.\",    \"favoriteFruit\": \"strawberry\"  },  {    \"_id\": \"5a0ae30d47367b43a9f4e819\",    \"index\": 3,    \"guid\": \"8192649f-8fc0-481d-a32a-50c6919d10c0\",    \"isActive\": true,    \"balance\": \"$3,221.56\",    \"picture\": \"http://placehold.it/32x32\",    \"age\": 28,    \"eyeColor\": \"blue\",    \"name\": \"Cash Cannon\",    \"gender\": \"male\",    \"company\": \"ACCUPRINT\",    \"email\": \"cashcannon@accuprint.com\",    \"phone\": \"+1 (990) 410-3156\",    \"address\": \"454 George Street, Austinburg, Iowa, 8682\",    \"about\": \"Consequat consectetur minim est in incididunt. Sunt enim voluptate amet consectetur voluptate enim id. Laborum incididunt pariatur consequat qui sint nisi aliquip aliquip consequat nostrud excepteur. Commodo ipsum officia sint pariatur voluptate officia minim occaecat eiusmod commodo ut ipsum consequat excepteur.\",    \"registered\": \"2015-07-11T02:59:53 -02:00\",    \"latitude\": 33.528547,    \"longitude\": 57.927556,    \"tags\": [      \"officia\",      \"labore\",      \"sit\",      \"mollit\",      \"dolore\",      \"labore\",      \"elit\"    ],    \"friends\": [      {        \"id\": 0,        \"name\": \"Burch Cherry\"      },      {        \"id\": 1,        \"name\": \"Bray Roberts\"      },      {        \"id\": 2,        \"name\": \"Bridgette Hester\"      }    ],    \"greeting\": \"Hello, Cash Cannon! You have 3 unread messages.\",    \"favoriteFruit\": \"strawberry\"  },  {    \"_id\": \"5a0ae30d97744a9d693c381b\",    \"index\": 4,    \"guid\": \"71e9af8c-8c30-4810-b4d0-92a4c6fd434c\",    \"isActive\": true,    \"balance\": \"$3,717.56\",    \"picture\": \"http://placehold.it/32x32\",    \"age\": 29,    \"eyeColor\": \"brown\",    \"name\": \"Jenny Mendoza\",    \"gender\": \"female\",    \"company\": \"DUFLEX\",    \"email\": \"jennymendoza@duflex.com\",    \"phone\": \"+1 (864) 549-3766\",    \"address\": \"578 Autumn Avenue, Stagecoach, California, 1739\",    \"about\": \"Ipsum ut aliqua qui nisi quis incididunt dolore fugiat. Lorem qui pariatur occaecat reprehenderit laboris nulla eiusmod. Ut sunt sint nostrud incididunt anim sint. Aliquip magna velit enim exercitation proident qui enim amet ullamco pariatur commodo Lorem. Veniam amet in aliqua mollit sint sit ex.\",    \"registered\": \"2014-09-25T09:34:35 -02:00\",    \"latitude\": -1.133185,    \"longitude\": -77.004934,    \"tags\": [      \"laborum\",      \"officia\",      \"eu\",      \"aute\",      \"et\",      \"aliqua\",      \"eu\"    ],    \"friends\": [      {        \"id\": 0,        \"name\": \"Marla Macias\"      },      {        \"id\": 1,        \"name\": \"Lauri Adkins\"      },      {        \"id\": 2,        \"name\": \"Coleman Keller\"      }    ],    \"greeting\": \"Hello, Jenny Mendoza! You have 1 unread messages.\",    \"favoriteFruit\": \"banana\"  }]";
  jsonStringToGenericContainer( a_json, gc );
  cout << "The container from the json string is:\n";
  gc.dump( cout );
  cout << "\n\nAll done Folks!!!\n\n";
  return 0;
}
