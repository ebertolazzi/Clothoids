/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2018                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Paolo Bevilacqua and Enrico Bertolazzi                              |
 |                                                                          |
 |      (1) Dipartimento di Ingegneria e Scienza dell'Informazione          |
 |      (2) Dipartimento di Ingegneria Industriale                          |
 |                                                                          |
 |      Universita` degli Studi di Trento                                   |
 |      email: paolo.bevilacqua@unitn.it                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: BBox.hxx
///

namespace G2lib {

  using std::setw;
  using std::vector;
  using std::pair;

  using std::make_shared;
  using std::shared_ptr; // promemoria shared_ptr<Foo>(&foo, [](void*){});

  /*\
   |   ____  ____
   |  | __ )| __ )  _____  __
   |  |  _ \|  _ \ / _ \ \/ /
   |  | |_) | |_) | (_) >  <
   |  |____/|____/ \___/_/\_\
  \*/
  //!
  //! Class to manipulate bounding box
  //!
  class BBox {
  public:
    typedef shared_ptr<BBox const> PtrBBox;
    typedef pair<PtrBBox,PtrBBox>  PairPtrBBox;
    typedef vector<PtrBBox>        VecPtrBBox;
    typedef vector<PairPtrBBox>    VecPairPtrBBox;

  private:
    real_type m_xmin; //!< left bottom
    real_type m_ymin; //!< left bottom
    real_type m_xmax; //!< right top
    real_type m_ymax; //!< right top
    int_type  m_id;   //!< id of the bbox
    int_type  m_ipos; //!< rank of the bounding box used in external algorithms

    BBox();

    BBox( BBox const & ) = default;
    BBox( BBox && ) = default;

  public:

    //!
    //! Construct a bounding box with additional information
    //!
    //! \param[in] xmin x-minimimum box coordinate
    //! \param[in] ymin y-minimimum box coordinate
    //! \param[in] xmax x-maximum box coordinate
    //! \param[in] ymax y-maximum box coordinate
    //! \param[in] id   identifier of the box
    //! \param[in] ipos ranking position of the box
    //!
    BBox(
      real_type xmin,
      real_type ymin,
      real_type xmax,
      real_type ymax,
      int_type  id,
      int_type  ipos
    ) {
      m_xmin = xmin;
      m_ymin = ymin;
      m_xmax = xmax;
      m_ymax = ymax;
      m_id   = id;
      m_ipos = ipos;
    }

    //!
    //! Build a buonding box that cover a list of bounding box
    //!
    //! \param[in] bboxes list of bounding box
    //! \param[in] id     identifier of the box
    //! \param[in] ipos   ranking position of the box
    //!
    BBox(
      vector<PtrBBox> const & bboxes,
      int_type                id,
      int_type                ipos
    ) {
      m_id   = id;
      m_ipos = ipos;
      this -> join( bboxes );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    real_type Xmin() const { return m_xmin; } //!< x-minimum coordinate of the bbox
    real_type Ymin() const { return m_ymin; } //!< y-minimum coordinate of the bbox
    real_type Xmax() const { return m_xmax; } //!< x-maximum coordinate of the bbox
    real_type Ymax() const { return m_ymax; } //!< y-maximum coordinate of the bbox

    int_type const & Id()   const { return m_id; }   //!< return BBOX id
    int_type const & Ipos() const { return m_ipos; } //!< return BBOX position

    //!
    //! copy a bbox
    //!
    BBox const &
    operator = ( BBox const & rhs ) {
      m_xmin = rhs.m_xmin;
      m_ymin = rhs.m_ymin;
      m_xmax = rhs.m_xmax;
      m_ymax = rhs.m_ymax;
      m_id   = rhs.m_id;
      m_ipos = rhs.m_ipos;
      return *this;
    }

    //!
    //! detect if two bbox collide
    //!
    bool
    collision( BBox const & box ) const {
      return !( (box.m_xmin > m_xmax ) ||
                (box.m_xmax < m_xmin ) ||
                (box.m_ymin > m_ymax ) ||
                (box.m_ymax < m_ymin ) );
    }

    //!
    //! Build bbox for a list of bbox
    //!
    void
    join( vector<PtrBBox> const & bboxes );

    //!
    //! distance of the point `(x,y)` to the bbox
    //!
    real_type
    distance( real_type x, real_type y ) const;

    //!
    //! Maximum distance of the point `(x,y)` to the point of bbox
    //!
    real_type
    maxDistance( real_type x, real_type y ) const;

    //!
    //! Pretty print a bbox
    //!
    void
    print( ostream_type & stream ) const {
      fmt::print( stream,
        "BBOX (xmin,ymin,xmax,ymax) = ( {}, {}, {}, {} )\n",
        m_xmin, m_ymin, m_xmax, m_ymax
      );
    }
  };

  //!
  //! Pretty print a bbox
  //!
  inline
  ostream_type &
  operator << ( ostream_type & stream, BBox const & bb ) {
    bb.print(stream);
    return stream;
  }

}

///
/// eof: BBox.hxx
///
