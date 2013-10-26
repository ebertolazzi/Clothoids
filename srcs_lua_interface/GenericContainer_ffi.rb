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
  
  HOST_OS = RbConfig::CONFIG['host_os']
  
  case HOST_OS
  when /mac|darwin/
    ffi_lib %w[libGenericContainer.dylib ./libGenericContainer.dylib ./libs/libGenericContainer.dylib]
  when /linux|cygwin|bsd/
    ffi_lib %w[libGenericContainer.so ./libGenericContainer.so ./libs/libGenericContainer.so]
  when /mswin|win|mingw/
    ffi_lib %w[libGenericContainer.dll ./libGenericContainer.dll ./libs/libGenericContainer.dll]
  else
  end

  attach_function :GC_new,                        [ :string ], :int
  attach_function :GC_select,                     [ :string ], :int
  attach_function :GC_delete,                     [ :string ], :int
  attach_function :GC_get_type,                   [], :int
  attach_function :GC_get_type_name,              [], :string
  attach_function :GC_print,                      [], :int
  attach_function :GC_mem_ptr,                    [ :string ], :pointer
  # head movement
  attach_function :GC_pop_head,                   [], :int
  attach_function :GC_reset_head,                 [], :int
  # bool
  attach_function :GC_set_bool,                   [ :int ], :int
  attach_function :GC_get_bool,                   [], :int
  attach_function :GC_set_vector_of_bool,         [ :pointer, :int ], :int
  attach_function :GC_set_empty_vector_of_bool,   [], :int
  attach_function :GC_push_bool,                  [ :int ], :int
  attach_function :GC_get_bool_at_pos,            [ :int ], :int
  # int
  attach_function :GC_set_int,                    [ :int ], :int
  attach_function :GC_get_int,                    [], :int
  attach_function :GC_set_vector_of_int,          [ :pointer, :int ], :int
  attach_function :GC_set_empty_vector_of_int,    [], :int
  attach_function :GC_push_int,                   [ :int ], :int
  attach_function :GC_get_int_at_pos,             [ :int ], :int
  # real
  attach_function :GC_set_real,                   [ :double ], :int
  attach_function :GC_get_real,                   [], :double
  attach_function :GC_set_vector_of_real,         [ :pointer, :int ], :int
  attach_function :GC_set_empty_vector_of_real,   [], :int
  attach_function :GC_push_real,                  [ :double ], :int
  attach_function :GC_get_real_at_pos,            [ :int ], :double
  # string
  attach_function :GC_set_string,                 [ :string ], :int
  attach_function :GC_get_string,                 [], :string
  attach_function :GC_set_vector_of_string,       [ :pointer, :int ], :int
  attach_function :GC_set_empty_vector_of_string, [], :int
  attach_function :GC_push_string,                [ :string ], :int
  attach_function :GC_get_string_at_pos,          [ :int ], :string
  # generic vector
  attach_function :GC_set_vector,                 [ :int ], :int
  attach_function :GC_set_empty_vector,           [], :int
  attach_function :GC_set_vector_position,        [ :int ], :int
  attach_function :GC_get_vector_size,            [], :int
  # generic map
  attach_function :GC_set_map,                    [], :int
  attach_function :GC_get_map,                    [], :int
  attach_function :GC_get_key,                    [], :string
  attach_function :GC_set_map_position,           [ :string ], :int
  
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
      self.GC_set_vector_position i
      res << self.GC_to_hash
      self.GC_pop_head
    end
    return res
  end

  def self.toMixedHash
    self.GC_get_map
    keys = []
    while (key = self.GC_get_key) do keys << key end
    res = {}
    keys.each do |key|
      self.GC_set_map_position key
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
  def self.fromBoolVector(var)
    if var.all? { |v| (v.class == :TrueClass) || (v.class == :FalseClass) } then
      self.GC_set_empty_vector_of_bool
      var.each { |v| self.GC_push_bool v }
      self.GC_pop_head
    else
      self.fromIntVector var
    end
  end

  def self.fromIntVector(var)
    if var.all? { |v| v.class == Fixnum } then
      self.GC_set_empty_vector_of_int
      var.each { |v| self.GC_push_int v }
    else
      self.fromRealVector var
    end
  end

  def self.fromRealVector(var)
    if var.all? { |v| v.class == Float } then
      self.GC_set_empty_vector_of_real
      var.each { |v| self.GC_push_real v }
    else
      self.fromStringVector var
    end
  end

  def self.fromStringVector(var)
    if var.all? { |v| v.class == String } then
      self.GC_set_empty_vector_of_string
      var.each { |v| self.GC_push_string v }
    else
      self.fromVector var
    end
  end

  def self.fromVector(var)
    nelem = var.length
    self.GC_set_vector nelem
    (0..nelem-1).each do |i|
      self.GC_set_vector_position i
      self.GC_from_hash var[i]
      self.GC_pop_head
    end
  end

  def self.fromMixedHash(var)
    self.GC_set_map
    var.each do |k,v|
      self.GC_set_map_position k.to_s
      self.GC_from_hash v
      self.GC_pop_head
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
    elsif var.kind_of?(Array) then
      tmp = var[0]
      if tmp.kind_of?(TrueClass) || tmp.kind_of?(FalseClass) then
        self.fromBoolVector var
      elsif tmp.kind_of?(Fixnum) then
        self.fromIntVector var
      elsif tmp.kind_of?(Float) then
        self.fromRealVector var
      elsif tmp.kind_of?(String) then
        self.fromStringVector var
      else
        self.fromVector var
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
      ::GenericContainer.GC_new @id
      ObjectSpace.define_finalizer(self, proc { ::GenericContainer.GC_delete @id })
      return @id
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
