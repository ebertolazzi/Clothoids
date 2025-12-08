#ifdef UTILS_OS_LINUX

#include "Utils_fmt.hh"

class mystream : public std::streambuf
{
protected:
  virtual std::streamsize
  xsputn( const char * s, std::streamsize n ) override
  {
    mexPrintf( "%.*s", n, s );
    mexEvalString( "drawnow;" );
    return n;
  }

  virtual int
  overflow( int c = EOF ) override
  {
    if ( c != EOF )
    {
      mexPrintf( "%.1s", &c );
      mexEvalString( "drawnow;" );
    }
    return 1;
  }
};

class scoped_redirect_cout
{
public:
  scoped_redirect_cout()
  {
    old_buf = std::cout.rdbuf();
    std::cout.rdbuf( &mout );

    // Crea un file stream che punti a MATLAB
    matlab_stream = fopen( "matlab_output", "w" );
    old_stdout    = stdout;
    stdout        = matlab_stream;
  }

  ~scoped_redirect_cout()
  {
    std::cout.rdbuf( old_buf );
    if ( matlab_stream )
    {
      fclose( matlab_stream );
      stdout = old_stdout;
    }
  }

private:
  mystream         mout;
  std::streambuf * old_buf;
  FILE *           matlab_stream;
  FILE *           old_stdout;
};

// Funzione personalizzata per l'output di fmt
void
matlab_output_handler( const char * buffer, size_t size )
{
  mexPrintf( "%.*s", static_cast<int>( size ), buffer );
  mexEvalString( "drawnow;" );
}

// Specializzazione per utilizzare il nostro handler
template <>
struct fmt::formatter<matlab_output_handler>
{
  constexpr auto
  parse( format_parse_context & ctx )
  {
    return ctx.begin();
  }

  template <typename FormatContext>
  auto
  format( const matlab_output_handler &, FormatContext & ctx )
  {
    return ctx.out();
  }
};

static scoped_redirect_cout mycout_redirect;

#endif
