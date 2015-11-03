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

require 'ffi'

module GenericContainer
  
  extend FFI::Library

  @GENERIC_CONTAINER_OK        = 0 ;
  @GENERIC_CONTAINER_BAD_TYPE  = 1 ;
  @GENERIC_CONTAINER_NO_DATA   = 2 ;
  @GENERIC_CONTAINER_NOT_EMPTY = 3 ;
  @GENERIC_CONTAINER_BAD_HEAD  = 4 ;

  @error_messages = [ "Generic Container, ok",
                      "Generic Container, bad type",
                      "Generic Container, no data",
                      "Generic Container, not empty",
                      "Generic Container, bad head" ] ;

  if @libGenericContainer then
    ffi_lib @libGenericContainer
  else
    HOST_OS  = RbConfig::CONFIG['host_os']
    @libname = "GenericContainer"
    @ext     = ".noextension"
    case HOST_OS
    when /mac|darwin/
      @ext = ".dylib" ;
    when /linux|cygwin|bsd/
      @ext = ".so" ;
    when /mswin|win|mingw/
      @ext = ".dll" ;
    end
    ffi_lib [ @libname+@ext,
              "./lib/"+@libname+@ext,
              "../lib/"+@libname+@ext,
              "./libs/"+@libname+@ext,
              "../libs/"+@libname+@ext,
              "lib"+@libname+@ext,
              "./lib/lib"+@libname+@ext,
              "../lib/lib"+@libname+@ext,
              "./libs/lib"+@libname+@ext,
              "../libs/lib"+@libname+@ext ]
  end

  attach_function :GC_select,                     [ :string ], :int
  attach_function :GC_delete,                     [ :string ], :int
  attach_function :GC_fill_for_test,              [ :string ], :int

  # head movement
  attach_function :GC_pop_head,                   [], :int
  attach_function :GC_reset_head,                 [], :int

  attach_function :GC_get_type,                   [], :int
  attach_function :GC_get_type_name,              [], :string
  attach_function :GC_print,                      [], :int
  attach_function :GC_mem_ptr,                    [ :string ], :pointer

  # set
  attach_function :GC_set_bool,                   [ :int ], :int
  attach_function :GC_set_int,                    [ :int ], :int
  attach_function :GC_set_real,                   [ :double ], :int
  attach_function :GC_set_complex2,               [ :double, :double ], :int
  attach_function :GC_set_string,                 [ :string ], :int

  # get
  attach_function :GC_get_bool,                   [], :int
  attach_function :GC_get_int,                    [], :int
  attach_function :GC_get_long,                   [], :long
  attach_function :GC_get_real,                   [], :double
  attach_function :GC_get_complex_re,             [], :double
  attach_function :GC_get_complex_im,             [], :double
  attach_function :GC_get_string,                 [], :string
  
  # push
  attach_function :GC_push_bool,                  [ :int ], :int
  attach_function :GC_push_int,                   [ :int ], :int
  attach_function :GC_push_real,                  [ :double ], :int
  attach_function :GC_push_complex2,              [ :double, :double ], :int
  attach_function :GC_push_string,                [ :string ], :int

  # get with position
  attach_function :GC_get_bool_at_pos,            [ :int ], :int
  attach_function :GC_get_int_at_pos,             [ :int ], :int
  attach_function :GC_get_real_at_pos,            [ :int ], :double
  attach_function :GC_get_complex_real_at_pos,    [ :int ], :double
  attach_function :GC_get_complex_imag_at_pos,    [ :int ], :double
  attach_function :GC_get_string_at_pos,          [ :int ], :string

  attach_function :GC_get_matrix_num_rows,        [], :int
  attach_function :GC_get_matrix_num_cols,        [], :int

  attach_function :GC_get_real_at_coor,           [ :int, :int ], :double
  attach_function :GC_get_complex_real_at_coor,   [ :int, :int ], :double
  attach_function :GC_get_complex_imag_at_coor,   [ :int, :int ], :double

  # empty vectors
  attach_function :GC_set_empty_vector_of_bool,   [], :int
  attach_function :GC_set_empty_vector_of_int,    [], :int
  attach_function :GC_set_empty_vector_of_real,   [], :int
  attach_function :GC_set_empty_vector_of_complex,[], :int
  attach_function :GC_set_empty_vector_of_string, [], :int

  # vectors not used
  #attach_function :GC_set_vector_of_bool,         [ :pointer, :int ], :int
  #attach_function :GC_set_vector_of_int,          [ :pointer, :int ], :int
  #attach_function :GC_set_vector_of_real,         [ :pointer, :int ], :int
  #attach_function :GC_set_vector_of_complex,      [ :pointer, :int ], :int
  #attach_function :GC_set_vector_of_string,       [ :pointer, :int ], :int

  # generic vector
  attach_function :GC_set_vector,                 [ :int ], :int
  attach_function :GC_set_empty_vector,           [], :int
  attach_function :GC_push_vector_position,       [ :int ], :int
  attach_function :GC_get_vector_size,            [], :int

  # generic map
  attach_function :GC_set_map,                    [], :int
  attach_function :GC_init_map_key,               [], :int
  attach_function :GC_get_next_key,               [], :string
  attach_function :GC_push_map_position,          [ :string ], :int
  
  #
  #   ____  ____    _           _               _     
  #  / ___|/ ___|  | |_ ___    | |__   __ _ ___| |__  
  # | |  _| |      | __/ _ \   | '_ \ / _` / __| '_ \ 
  # | |_| | |___   | || (_) |  | | | | (_| \__ \ | | |
  #  \____|\____|___\__\___/___|_| |_|\__,_|___/_| |_|
  #            |_____|    |_____|
  #
  def self.toBoolVector
    nelem = self.GC_get_vector_size
    res = []
    (0..nelem-1).each { |i| res << self.GC_get_bool_at_pos(i) }
    return res
  end

  def self.toIntVector
    nelem = self.GC_get_vector_size
    res = []
    (0..nelem-1).each { |i| res << self.GC_get_int_at_pos(i) }
    return res
  end

  def self.toRealVector
    nelem = self.GC_get_vector_size
    res = []
    (0..nelem-1).each { |i| res << self.GC_get_real_at_pos(i) }
    return res
  end

  def self.toComplexVector
    nelem = self.GC_get_vector_size
    res = []
    (0..nelem-1).each { |i| res << [ self.GC_get_complex_real_at_pos(i), self.GC_get_complex_imag_at_pos(i) ] }
    return res
  end

  def self.toStringVector
    nelem = self.GC_get_vector_size
    res = []
    (0..nelem-1).each { |i| res << self.GC_get_string_at_pos(i) }
    return res
  end

  def self.toVector
    nelem = self.GC_get_vector_size
    res = []
    (0..nelem-1).each do |i|
      ok = self.GC_push_vector_position i
      raise RuntimeError, @error_messages[ok] unless ok == 0
      res << self.GC_to_hash
      ok = self.GC_pop_head
      raise RuntimeError, @error_messages[ok] unless ok == 0
    end
    return res
  end

  def self.toRealMatrix
    nr  = self.GC_get_matrix_num_rows
    nc  = self.GC_get_matrix_num_cols
    res = []
    (0...nr).each do |i|
      row = []
      (0...nc).each do |j|
        row << self.GC_get_real_at_coor(i,j) 
      end
      res << row
    end
    return res
  end

  def self.toComplexMatrix
    nr  = self.GC_get_matrix_num_rows
    nc  = self.GC_get_matrix_num_cols
    res = []
    (0...nr).each do |i|
      row = []
      (0...nc).each do |j|
        row << [ self.GC_get_complex_real_at_coor(i,j), self.GC_get_complex_imag_at_coor(i,j) ]
      end
      res << row
    end
    return res
  end

  def self.toMixedHash
    #self.GC_get_map
    keys = []
    self.GC_init_map_key
    while (key = self.GC_get_next_key) do keys << key end
    res = {}
    keys.each do |key|
      self.GC_push_map_position key
      res[key.to_sym] = self.GC_to_hash
      self.GC_pop_head
    end
    return res
  end

  def self.GC_to_hash
    type = self.GC_get_type_name
    case type
    #when /NOTYPE/
    #when /pointer/
    when /^bool_type$/
      res = self.GC_get_bool
    when /^int_type$/
      res = self.GC_get_int
    when /^long_type$/
      res = self.GC_get_long
    when /^real_type$/
      res = self.GC_get_real
    when /^string_type$/
      res = self.GC_get_string
    when /^vec_bool_type$/
      res = self.toBoolVector
    when /^vec_int_type$/
      res = self.toIntVector
    when /^vec_real_type$/
      res = self.toRealVector
    when /^vec_string_type$/
      res = self.toStringVector
    when /^vec_complex_type$/
      res = self.toComplexVector
    when /^mat_real_type$/
      res = self.toRealMatrix
    when /^mat_complex_type$/
      res = self.toComplexMatrix
    when /^vector_type$/
      res = self.toVector
    when /^map_type$/
      res = self.toMixedHash
    else
      puts "GC_to_hash #{self.GC_get_type_name} not managed"
      res = nil
    end
    return res
  end
  #
  #    ____  ____      __                         _               _     
  #   / ___|/ ___|    / _|_ __ ___  _ __ ___     | |__   __ _ ___| |__  
  #  | |  _| |       | |_| '__/ _ \| '_ ` _ \    | '_ \ / _` / __| '_ \ 
  #  | |_| | |___    |  _| | | (_) | | | | | |   | | | | (_| \__ \ | | |
  #   \____|\____|___|_| |_|  \___/|_| |_| |_|___|_| |_|\__,_|___/_| |_|
  #             |_____|                     |_____|                     
  #
  def self.addToVector(vec)
    # empty vectors
    if vec[0].kind_of?(TrueClass) or vec[0].kind_of?(FalseClass) then
      self.GC_set_empty_vector_of_bool
    elsif vec[0].kind_of?(Fixnum) then
      self.GC_set_empty_vector_of_int
    elsif vec[0].kind_of?(Float) then
      self.GC_set_empty_vector_of_real
    elsif vec[0].kind_of?(String) then
      self.GC_set_empty_vector_of_string
    elsif vec[0].kind_of?(Symbol) then
      self.GC_set_empty_vector_of_string
    else
      self.fromVector vec
      return ;
    end
    # add elements
    vec.each { |var|
      if var.kind_of?(TrueClass) then
        self.GC_push_bool 1
      elsif var.kind_of?(FalseClass) then
        self.GC_push_bool 0
      elsif var.kind_of?(Fixnum) then
        self.GC_push_int var
      elsif var.kind_of?(Float) then
        self.GC_push_real var
      elsif var.kind_of?(String) then
        self.GC_push_string var
      elsif var.kind_of?(Symbol) then
        self.GC_push_string var.to_s
      else
        self.fromVector vec
        return ;
      end
    }
  end

  def self.fromVector(var)
    nelem = var.length
    self.GC_set_vector nelem
    (0..nelem-1).each do |i|
      ok = self.GC_push_vector_position i
      raise RuntimeError, @error_messages[ok] unless ok == 0
      self.GC_from_hash var[i]
      ok = self.GC_pop_head
      raise RuntimeError, @error_messages[ok] unless ok == 0
    end
  end

  def self.fromMixedHash(var)
    self.GC_set_map
    var.each do |k,v|
      ok = self.GC_push_map_position k.to_s
      raise RuntimeError, @error_messages[ok] unless ok == 0
      self.GC_from_hash v
      ok = self.GC_pop_head
      raise RuntimeError, @error_messages[ok] unless ok == 0
    end
  end

  def self.GC_from_hash(var)
    if var.kind_of?(TrueClass) then
      self.GC_set_bool 1
    elsif var.kind_of?(FalseClass) then
      self.GC_set_bool 0
    elsif var.kind_of?(Fixnum) then
      self.GC_set_int var
    elsif var.kind_of?(Float) then
      self.GC_set_real var
    elsif var.kind_of?(String) then
      self.GC_set_string var
    elsif var.kind_of?(Symbol) then
      self.GC_set_string var.to_s
    elsif var.kind_of?(Array) then
      if var.length > 0 then
        self.addToVector(var)
      else
        self.GC_set_empty_vector_of_bool
      end
    elsif var.kind_of?(Hash) then
      self.fromMixedHash var
    elsif var.respond_to? :to_hash then
      self.fromMixedHash var.to_hash
    elsif var.respond_to? :[] then
      self.fromMixedHash var
    else
      puts "Found MRB_TT_OBJECT: " ;
    end
  end

  class GenericContainer
    attr_reader :id
 
    def initialize
      @id = self.__id__.to_s
      ::GenericContainer.GC_select @id
      ObjectSpace.define_finalizer(self, proc { ::GenericContainer.GC_delete @id })
      return @id
    end

    def fill_for_test
      ::GenericContainer.GC_fill_for_test @id
    end

    def load(data)
      ::GenericContainer.GC_select    @id
      ::GenericContainer.GC_from_hash data
    end
    
    def get_pointer
      return ::GenericContainer.GC_mem_ptr @id
    end

    def get_data
      ::GenericContainer.GC_select @id
      return ::GenericContainer.GC_to_hash
    end

    def print
      ::GenericContainer.GC_select @id
      ::GenericContainer.GC_print
    end

  end

end
