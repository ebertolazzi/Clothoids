Introduction
============

STL iostream implementation using the library zlib.
This means that you can easily manipulate zipped
streams like any other STL ostream/istream.

To give you an idea, consider following snippet that prints "Hello World":

    ostringstream output_buffer;
    // writing data
    output_buffer << "Hello world" << endl ;

Now, the same snippet but with zipped output using zlib:

    // zip_ostream uses output_buffer as output buffer :)
    ozstream zipper( output_buffer );
    
    // writing data as usual
    zipper << "Hello world" << endl;

Or, to create gzipped files:

    ofstream file("hello_world.txt.gz");
    ogzstream gzfile(file);
    gzfile << "Hello world " << endl;

As you can see adding zipped buffers into your existing 
applications is quite straightforward.
To summarize, let's see some quick facts about zstream:

    * STL compliant,
    * any-stream-to-any-stream support,
    * char, wchar_t support,
    * fining tuning of compression properties,
    * support custom allocators (New!)


Based on the work of Jonathan de Halleux, published on 
CodeProject http://www.codeproject.com/Articles/4457/zipstream-bzip-stream-iostream-wrappers-for-the-zl