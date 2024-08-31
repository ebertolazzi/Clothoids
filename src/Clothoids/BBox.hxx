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
    using PtrBBox        = shared_ptr<BBox const>;
    using PairPtrBBox    = pair<PtrBBox,PtrBBox>;
    using VecPtrBBox     = vector<PtrBBox>;
    using VecPairPtrBBox = vector<PairPtrBBox>;

  private:
    real_type m_bbox[4]{0,0,0,0}; //!< [ xmin ymin xmax ymax ]
    integer   m_id{0};            //!< id of the bbox
    integer   m_ipos{0};          //!< rank of the bounding box used in external algorithms

    BBox();

    BBox( BBox const & ) = default;
    BBox( BBox && ) = default;

  public:

    //!
    //! Construct a bounding box with additional information
    //!
    //! \param[in] xmin \f$ x \f$-minimimum box coordinate
    //! \param[in] ymin \f$ y \f$-minimimum box coordinate
    //! \param[in] xmax \f$ x \f$-maximum box coordinate
    //! \param[in] ymax \f$ y \f$-maximum box coordinate
    //! \param[in] id   identifier of the box
    //! \param[in] ipos ranking position of the box
    //!
    BBox(
      real_type xmin,
      real_type ymin,
      real_type xmax,
      real_type ymax,
      integer   id,
      integer   ipos
    ) {
      m_bbox[0] = xmin;
      m_bbox[1] = ymin;
      m_bbox[2] = xmax;
      m_bbox[3] = ymax;
      m_id      = id;
      m_ipos    = ipos;
    }

    //!
    //! Construct a bounding box with additional information
    //!
    //! \param[in] bbox bbox [pmin, pmax]
    //! \param[in] id   identifier of the box
    //! \param[in] ipos ranking position of the box
    //!
    BBox(
      real_type const bbox[4],
      integer         id,
      integer         ipos
    ) {
      std::copy_n( bbox, 4, m_bbox );
      m_id   = id;
      m_ipos = ipos;
    }

    //!
    //! Construct a bounding box with additional information
    //!
    //! \param[in] bbox_min bounding box lower corner
    //! \param[in] bbox_max bounding box upper corner
    //! \param[in] id       identifier of the box
    //! \param[in] ipos     ranking position of the box
    //!
    BBox(
      real_type const bbox_min[2],
      real_type const bbox_max[2],
      integer         id,
      integer         ipos
    ) {
      std::copy_n( bbox_min, 2, m_bbox   );
      std::copy_n( bbox_max, 2, m_bbox+2 );
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
      integer                 id,
      integer                 ipos
    ) {
      m_id   = id;
      m_ipos = ipos;
      this -> join( bboxes );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    real_type const * bbox()     { return m_bbox; }
    real_type const * bbox_min() { return m_bbox; }
    real_type const * bbox_max() { return m_bbox+2; }

    real_type & x_min() { return m_bbox[0]; } //!< x-minimum coordinate of the bbox
    real_type & y_min() { return m_bbox[1]; } //!< y-minimum coordinate of the bbox
    real_type & x_max() { return m_bbox[2]; } //!< x-maximum coordinate of the bbox
    real_type & y_max() { return m_bbox[3]; } //!< y-maximum coordinate of the bbox

    real_type const & x_min() const { return m_bbox[0]; } //!< x-minimum coordinate of the bbox
    real_type const & y_min() const { return m_bbox[1]; } //!< y-minimum coordinate of the bbox
    real_type const & x_max() const { return m_bbox[2]; } //!< x-maximum coordinate of the bbox
    real_type const & y_max() const { return m_bbox[3]; } //!< y-maximum coordinate of the bbox

    integer const & Id()   const { return m_id; }   //!< return BBOX id
    integer const & Ipos() const { return m_ipos; } //!< return BBOX position

    //!
    //! copy a bbox
    //!
    BBox const &
    operator = ( BBox const & rhs ) {
      std::copy_n( rhs.m_bbox, 4, m_bbox );
      m_id   = rhs.m_id;
      m_ipos = rhs.m_ipos;
      return *this;
    }

    //!
    //! detect if two bbox collide
    //!
    bool
    collision( BBox const & box ) const {
      return !( (box.x_min() > this->x_max() ) ||
                (box.x_max() < this->x_min() ) ||
                (box.y_min() > this->y_max() ) ||
                (box.y_max() < this->y_min() ) );
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
    max_distance( real_type x, real_type y ) const;

    //!
    //! Pretty print a bbox
    //!
    void print( ostream_type & stream ) const;

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //!  Print on strem the `BBox` object
  //!
  //!  \param stream the output stream
  //!  \param bb     an instance of `BBox` object
  //!  \return the output stream
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
