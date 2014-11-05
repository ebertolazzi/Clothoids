#--------------------------------------------------------------------------#
#                                                                          |
#  Copyright (C) 2013                                                      |
#                                                                          |
#         , __                 , __                                        |
#        /|/  \               /|/  \                                       |
#         | __/ _   ,_         | __/ _   ,_                                |
#         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
#         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
#                       /|                   /|                            |
#                       \|                   \|                            |
#                                                                          |
#      E.Bertolazzi                                                        |
#      Dipartimento di Ingegneria Industriale                              |
#      Universita` degli Studi di Trento                                   |
#      email: enrico.bertolazzi@unitn.it                                   |
#                                                                          |
#--------------------------------------------------------------------------#

require 'pp'
require './GenericContainer_ffi.rb'

gc = ::GenericContainer::GenericContainer.new

data = {
  :infoLevel => 2,
  :bc => {
    :initial_x => true,
    :initial_v => true,
    :final_x   => false,
    :final_v   => false,
  },
  :Parameters => {
    :alpha             => 2,
    :beta              => 1,
    :gm                => 3,
    :uaControlMaxValue => 5,
    :uaControlMinValue => "pippo",
    :vettore  => [1,2,3,3.4,8,9,10],
    :vettore1 => [1,2,3,3,8,9,10],
    :vettore2 => [1,2,"222",3.4,[1,2],9,10],
    :vettore3 => [1.1,1.2,1],
  },
  :Controls => {
    :uaControl => {
      :type      => "aaaa",
      :epsi      => 0.1e-2,
      :tolerance => 0.5e-2,
    },
    :ubControl => {
      :type      => "bbb",
      :epsi      => 0.1e-2,
      :tolerance => 0.5e-2,
    },
  }
}

#gc.load data
#gc.print

gc.fill_for_test
gc.print
hsh = gc.get_data
pp hsh

