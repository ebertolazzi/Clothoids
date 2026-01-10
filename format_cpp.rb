#!/usr/bin/env ruby
# Script per formattare file C/C++ con clang-format

require 'fileutils'
require 'tempfile'
require 'json'

# Sostituisci la sezione CLANG_FORMAT_CONFIG con questa versione:
CLANG_FORMAT_CONFIG = {
  BasedOnStyle: "Google",
  IndentWidth: 2,
  UseTab: false,
  TabWidth: 2,
  ColumnLimit: 120,

  # --- FORZA UN PARAMETRO PER RIGA ---
  AlignAfterOpenBracket: "AlwaysBreak",
  BinPackParameters: false,
  BinPackArguments: false,
  AllowAllParametersOfDeclarationOnNextLine: false,
  AllowAllArgumentsOnNextLine: false,
  AlwaysBreakAfterReturnType: "None",
  AlwaysBreakAfterDefinitionReturnType: "None",  # MODIFICATO: da "All" a "None"
  ContinuationIndentWidth: 2,

  # --- ALLINEAMENTO VERTICALE (EFFETTO TABELLA) ---
  AlignConsecutiveDeclarations: true,
  AlignConsecutiveAssignments: true,
  AlignOperands: true,
  PointerAlignment: "Middle",
  DerivePointerAlignment: false,

  # --- PARENTESI E SPAZI ---
  BreakBeforeBraces: "Allman",
  SpacesInParentheses: true,
  SpacesInAngles: false,
  SpacesInSquareBrackets: false,
  SpaceBeforeParens: "ControlStatements",
  SpaceAfterCStyleCast: true,
  SpaceAfterTemplateKeyword: true,
  SpaceBeforeCpp11BracedList: false,
  Cpp11BracedListStyle: false,

  # --- COSTRUTTORI ---
  BreakConstructorInitializers: "BeforeComma",
  ConstructorInitializerAllOnOneLineOrOnePerLine: true,
  ConstructorInitializerIndentWidth: 2,

  # --- LOGICA DI INDENTAZIONE ---
  NamespaceIndentation: "All",
  FixNamespaceComments: true,
  AccessModifierOffset: -2,
  IndentCaseLabels: true,
  IndentPPDirectives: "None",
  
  # --- PULIZIA RIGHE ---
  MaxEmptyLinesToKeep: 2,
  KeepEmptyLinesAtTheStartOfBlocks: false,
  ReflowComments: true,
  SpacesBeforeTrailingComments: 2,
  SortIncludes: false,

  # --- FUNZIONI CORTE ---
  AllowShortFunctionsOnASingleLine: "InlineOnly",  # MODIFICATO: da "All" a "InlineOnly"
  AllowShortIfStatementsOnASingleLine: true,
  AllowShortLoopsOnASingleLine: true,
  AllowShortBlocksOnASingleLine: true,

  # --- PENALITÀ (FORZANO IL COMPORTAMENTO) ---
  # Alziamo la penalità per chi mette tutto su una riga
  PenaltyBreakAssignment: 100,
  PenaltyBreakBeforeFirstCallParameter: 1, 
  PenaltyExcessCharacter: 1000,
  PenaltyReturnTypeOnItsOwnLine: 1000,  # MODIFICATO: aumentata per evitare ritorni a capo

  # --- NUOVE IMPOSTAZIONI AGGIUNTE ---
  # Aggiungi queste impostazioni per gestire meglio le funzioni inline
  AllowShortLambdasOnASingleLine: "All",
  AllowShortCaseLabelsOnASingleLine: true,
  AlwaysBreakTemplateDeclarations: "No",
  
  # Per mantenere le funzioni inline su una riga
  BreakAfterJavaFieldAnnotations: false,
  BreakStringLiterals: false,
  
  # Aggiunto per gestire le funzioni nelle classi
  AllowShortFunctionsOnASingleLine: "Inline",  # MODIFICATO: imposta su "Inline"
  BraceWrapping: {
    AfterClass: false,
    AfterControlStatement: false,
    AfterEnum: false,
    AfterFunction: false,
    AfterNamespace: false,
    AfterObjCDeclaration: false,
    AfterStruct: false,
    AfterUnion: false,
    AfterExternBlock: false,
    BeforeCatch: false,
    BeforeElse: false,
    BeforeLambdaBody: false,
    BeforeWhile: false,
    IndentBraces: false,
    SplitEmptyFunction: true,
    SplitEmptyRecord: true,
    SplitEmptyNamespace: true
  },
  # Aggiungi queste se non sono già presenti:
  SeparateDefinitionBlocks: "Leave",
  EmptyLineAfterAccessModifier: "Never",
  CompactNamespaces: true,
}.freeze


# Estensioni dei file da processare
EXTENSIONS = ['.cc', '.hh', '.cxx', '.hxx']

# Directory da saltare (aggiunto)
EXCLUDED_DIRS = ['3rd','Eigen'].freeze

# Colori per output console
COLOR_RESET = "\e[0m"
COLOR_RED = "\e[31m"
COLOR_GREEN = "\e[32m"
COLOR_YELLOW = "\e[33m"
COLOR_BLUE = "\e[34m"

def log_info(message)
  puts "#{COLOR_BLUE}[INFO]#{COLOR_RESET} #{message}"
end

def log_success(message)
  puts "#{COLOR_GREEN}[SUCCESS]#{COLOR_RESET} #{message}"
end

def log_warning(message)
  puts "#{COLOR_YELLOW}[WARNING]#{COLOR_RESET} #{message}"
end

def log_error(message)
  puts "#{COLOR_RED}[ERROR]#{COLOR_RESET} #{message}"
end

def check_clang_format
  # Verifica se clang-format è installato
  unless system('which clang-format > /dev/null 2>&1')
    log_error "clang-format non è installato!"
    log_info "Per installarlo:"
    log_info "  Ubuntu/Debian: sudo apt-get install clang-format"
    log_info "  macOS: brew install clang-format"
    log_info "  Windows: scoop install llvm"
    exit 1
  end
  
  # Verifica la versione di clang-format
  version_output = `clang-format --version 2>&1`
  log_info "Trovato: #{version_output.split("\n").first}"
end

def create_temp_config_file
  yaml_file_path = File.join(Dir.tmpdir, "clang-format-config-#{Process.pid}.yml")
  
  # Converti la configurazione in formato YAML
  yaml_content = []
  
  CLANG_FORMAT_CONFIG.each do |key, value|
    case value
    when String
      yaml_content << "#{key}: #{value}"
    when TrueClass, FalseClass
      yaml_content << "#{key}: #{value.to_s.downcase}"
    when Array
      # Gestione speciale per IncludeCategories
      if key == :IncludeCategories
        yaml_content << "#{key}:"
        value.each do |category|
          yaml_content << "  - Regex: '#{category[:Regex]}'"
          yaml_content << "    Priority: #{category[:Priority]}"
          if category.key?(:SortPriority)
            yaml_content << "    SortPriority: #{category[:SortPriority]}"
          end
        end
      else
        yaml_content << "#{key}: #{value}"
      end
    when Hash
      # Gestione per hash (es. BraceWrapping se presente)
      yaml_content << "#{key}:"
      value.each do |subkey, subvalue|
        if subvalue.is_a?(TrueClass) || subvalue.is_a?(FalseClass)
          yaml_content << "  #{subkey}: #{subvalue.to_s.downcase}"
        else
          yaml_content << "  #{subkey}: #{subvalue}"
        end
      end
    else
      yaml_content << "#{key}: #{value}"
    end
  end
  
  File.write(yaml_file_path, yaml_content.join("\n"))
  yaml_file_path
end

def find_source_files(directory = '.')
  files = []
  
  EXTENSIONS.each do |ext|
    # Cerca file con l'estensione specificata (case insensitive)
    pattern = File.join(directory, '**', "*#{ext}")
    found = Dir.glob(pattern, File::FNM_CASEFOLD)
    
    # Filtra le directory escluse (MODIFICATO)
    found.reject! do |file|
      EXCLUDED_DIRS.any? { |excluded_dir| file.include?("/#{excluded_dir}/") }
    end
    
    files += found
    
    # Cerca anche versioni con lettere maiuscole
    pattern_upper = File.join(directory, '**', "*#{ext.upcase}")
    found_upper = Dir.glob(pattern_upper)
    
    # Filtra anche per le versioni maiuscole (MODIFICATO)
    found_upper.reject! do |file|
      EXCLUDED_DIRS.any? { |excluded_dir| file.include?("/#{excluded_dir}/") }
    end
    
    files += found_upper
  end
  
  # Rimuovi duplicati (nel caso siano trovati sia con case sensitive che insensitive)
  files.uniq.sort
end

def format_file(file_path, config_path)
  # Esegue clang-format sul file
  command = "clang-format -i -style=file:#{config_path} \"#{file_path}\""
  
  if system(command)
    log_success "Formattato: #{file_path}"
    return true
  else
    log_error "Errore nella formattazione di: #{file_path}"
    return false
  end
end

def dry_run(file_path, config_path)
  # Mostra le differenze senza applicare le modifiche
  temp_file = Tempfile.new(['clang-format', File.extname(file_path)])
  temp_file.close
  
  command = "clang-format -style=file:#{config_path} \"#{file_path}\" > \"#{temp_file.path}\""
  
  if system(command)
    diff = `diff -u "#{file_path}" "#{temp_file.path}" 2>&1`
    
    if diff.empty?
      log_info "Nessuna modifica necessaria per: #{file_path}"
    else
      log_warning "Modifiche necessarie per: #{file_path}"
      puts diff if ARGV.include?('--show-diff')
    end
    
    temp_file.unlink
    return !diff.empty?
  else
    log_error "Errore nella verifica di: #{file_path}"
    return false
  end
end

def print_summary(total, formatted, errors)
  puts "\n" + "=" * 50
  puts "RIEPILOGO:"
  puts "  File totali trovati: #{total}"
  puts "  File formattati: #{formatted}"
  puts "  Errori: #{errors}"
  puts "=" * 50
end

def print_usage
  puts "Uso: ruby #{File.basename(__FILE__)} [OPZIONI] [DIRECTORY]"
  puts ""
  puts "Opzioni:"
  puts "  --dry-run        Mostra cosa verrebbe modificato senza applicare cambiamenti"
  puts "  --show-diff      Mostra le differenze durante il dry-run"
  puts "  --help           Mostra questo messaggio di aiuto"
  puts ""
  puts "Directory:"
  puts "  Specifica la directory da scansionare (default: corrente)"
  puts ""
  puts "Estensioni processate:"
  EXTENSIONS.each { |ext| puts "  *#{ext}" }
  puts ""
  puts "Directory escluse:"
  EXCLUDED_DIRS.each { |dir| puts "  */#{dir}/" }
end

def main
  # Parsing degli argomenti
  dry_run_mode = false
  show_diff = false
  target_dir = '.'
  
  ARGV.each do |arg|
    case arg
    when '--dry-run'
      dry_run_mode = true
    when '--show-diff'
      show_diff = true
    when '--help', '-h'
      print_usage
      exit 0
    else
      # Verifica se è una directory valida
      if File.directory?(arg)
        target_dir = arg
      elsif arg.start_with?('-')
        log_error "Opzione sconosciuta: #{arg}"
        print_usage
        exit 1
      end
    end
  end
  
  log_info "Ricerca file in: #{File.expand_path(target_dir)}"
  log_info "Directory escluse: #{EXCLUDED_DIRS.join(', ')}"
  
  # Verifica clang-format
  check_clang_format
  
  # Crea file di configurazione temporaneo
  config_path = create_temp_config_file
  log_info "Configurazione clang-format creata: #{config_path}"
  
  # Trova i file sorgente
  files = find_source_files(target_dir)
  
  if files.empty?
    log_warning "Nessun file trovato con le estensioni specificate"
    File.delete(config_path) if File.exist?(config_path)
    exit 0
  end
  
  log_info "Trovati #{files.length} file da processare"
  
  # Processa i file
  formatted_count = 0
  error_count = 0
  
  files.each do |file|
    begin
      if dry_run_mode
        needs_formatting = dry_run(file, config_path)
        formatted_count += 1 if needs_formatting
      else
        success = format_file(file, config_path)
        if success
          formatted_count += 1
        else
          error_count += 1
        end
      end
    rescue => e
      log_error "Errore processando #{file}: #{e.message}"
      error_count += 1
    end
  end
  
  # Pulisci il file di configurazione temporaneo
  File.delete(config_path) if File.exist?(config_path)
  
  # Stampa il riepilogo
  print_summary(files.length, formatted_count, error_count)
  
  exit(error_count > 0 ? 1 : 0)
end

if __FILE__ == $0
  main
end
