#include <OverlapAreas.hxx>
#include <mpi.h>

namespace Engine
{
    OverlapAreas::OverlapAreas( int id, int nprocs, int size, int over )
    {
        _pid = id;
        _dim = sqrt( nprocs );
        if ( _dim*_dim != nprocs )
        {
            if ( id == 0 )
                std::cerr << "ABORTING PROGRAM:\n\tPlease, use a number of processes equal to NP*NP\n";
            abort( );
        }
        if ( size%_dim != 0 )
        {
            if ( id == 0 )
                std::cerr << "ABORTING PROGRAM:\n\tSize of the raster [" << size << "] is not multiple of N [" << _dim << "]\n";
            abort( );
        }
        int lsize = size/_dim;
        if ( lsize%2 != 0 )
        {
            if ( id == 0 )
                std::cerr << "ABORTING PROGRAM:\n\tLocal size of the raster [" << lsize << "] is not multiple of 2\n";
            abort( );
        }

        Rectangle<int> global = Rectangle<int>( Size<int>( size, size ), Point2D<int>( 0, 0 ) );
        _pos = Point2D<int>( id%_dim, id/_dim );
        _area = Rectangle<int>( Size<int>( lsize, lsize ), _pos*lsize );
        _boun = global.intersection( _area + Size<int>( 2*over, 2*over) - Point2D<int>(over,over) );

        // LEFT BORDER
        if ( isLeft( ) )
        {
            _l_neig = id - 1;
            _l_area = Rectangle<int>( _area.left( ), _area.top( ), _area.left( )+over-1, _area.bottom( ) );
            _l_boun = _l_area - Point2D<int>( over, 0 );
            _left._n = _l_neig;
            _left._local = _l_area;
            _left._bound = _l_boun;
            _topleft._n = _left._n;
            _topleft._local = _left._local.reSize( Size<int>( lsize/2, over ) );
            _topleft._bound = _left._bound.reSize( Size<int>( lsize/2, over ) );
         }
        MPI_Barrier( MPI_COMM_WORLD );

        // RIGHT BORDER
        if ( isRight( ) )
        {
            _r_neig = id + 1;
            _r_area = Rectangle<int>( _area.right( ) - over + 1, _area.top( ), _area.right( ), _area.bottom( ) );
            _r_boun = _r_area + Point2D<int>( over, 0 );
            _right._n = _r_neig;
            _right._local = _r_area;
            _right._bound = _r_boun;

            _topright._n = _right._n;
            _topright._local = _right._local.reSize( Size<int>( lsize/2, over ) );
            _topright._bound = _right._bound.reSize( Size<int>( lsize/2, over ) );
        }
        MPI_Barrier( MPI_COMM_WORLD );
        
        // TOP BORDER
        if ( isTop( ) )
        {
            _t_neig = id - _dim;
            _t_area = Rectangle<int>( _boun.left( ), _area.top( ), _boun.right( ), _area.top( ) + over - 1 );
            _t_boun = _t_area - Point2D<int>( 0, over );
            _top._n = id - _dim;;
            _top._local = _t_area;
            _top._bound = _t_boun;
        }
        MPI_Barrier( MPI_COMM_WORLD );

        // BOTTOM BORDER
        if ( isBottom( ) )
        {
            _b_neig = id + _dim;
            _b_area = Rectangle<int>( _boun.left( ), _area.bottom( ) - over + 1, _boun.right( ), _area.bottom( ) );
            _b_boun = _b_area + Point2D<int>( 0, over );
            _bottom._n = id + _dim;;
            _bottom._local = _b_area;
            _bottom._bound = _b_boun;
        }

        if ( _pos._x%2 == _pos._y%2 )
        {
            if ( isLeft( ) ) push_back( _l_neig, _l_area, _l_boun );
            if ( isRight( ) ) push_back( _r_neig, _r_area, _r_boun );
            if ( isTop( ) ) push_back( _t_neig, _t_area, _t_boun );
            if ( isBottom( ) ) push_back( _b_neig, _b_area, _b_boun );
        }
        else
        {
            if ( isRight( ) ) push_back( _r_neig, _r_area, _r_boun );
            if ( isLeft( ) ) push_back( _l_neig, _l_area, _l_boun );
            if ( isBottom( ) ) push_back( _b_neig, _b_area, _b_boun );
            if ( isTop( ) ) push_back( _t_neig, _t_area, _t_boun );
        }

#ifdef TMP
        if ( _worldPos._x < dim-1 )
        {
            _nearby[ind] = _id + 1;
            _nb_area[ind] = Rectangle<int>( _ownedArea.right( )-_overlap+1, _ownedArea.top( ), _ownedArea.right( ), _ownedArea.bottom( ) );
            _nb_exta[ind] = Rectangle<int>( _ownedArea.right( )+1, _ownedArea.top( ), _ownedArea.right( )+_overlap, _ownedArea.bottom( ) );
        }
        else
            _nearby[ind] = -1;

        // TOP BORDER
        ind = _worldPos._x%2 == _worldPos._y%2 ? 2 : 3;
        if ( _worldPos._y > 0 )
        {
            _nearby[ind] = _id - dim;
            _nb_area[ind] = Rectangle<int>( _boundaries.left( ), _ownedArea.top( ), _boundaries.right( ), _ownedArea.top( )+_overlap-1 );
            _nb_exta[ind] = Rectangle<int>( _boundaries.left( ), _boundaries.top( ), _boundaries.right( ), _ownedArea.top( )-1 );
        }
        else
            _nearby[ind] = -1;

        // BOTTOM BORDER
        ind = _worldPos._x%2 == _worldPos._y%2 ? 3 : 2;
        if ( _worldPos._y < dim-1 )
        {
            _nearby[ind] = _id + dim;
            _nb_area[ind] = Rectangle<int>( _boundaries.left( ), _ownedArea.bottom( )-_overlap+1, _boundaries.right( ), _ownedArea.bottom( ) );
            _nb_exta[ind] = Rectangle<int>( _boundaries.left( ), _ownedArea.bottom( )+1, _boundaries.right( ), _boundaries.bottom( ) );
        }
        else
            _nearby[ind] = -1;
#endif

    }

    void OverlapAreas::push_back( int neig, Rectangle<int> area, Rectangle<int> boun )
    {
        _v_neig.push_back( neig );
        _v_area.push_back( area );
        _v_boun.push_back( boun );
    }

    void OverlapAreas::abort( )
    {
        MPI_Barrier( MPI_COMM_WORLD );
        MPI_Finalize( );
        exit(0);
    }

#ifdef RGT
    OverlapAreas::OverlapAreas( int pid, int nprocs, int overlap, int gsize ) : _pid( pid ), _nprocs( nprocs ), _overlap( overlap ), _gsize( gsize )
    {
        // position of world related to the complete set of computer nodes
        _dim = sqrt( _nprocs );
        if ( _dim*_dim != _nprocs )
        {
            if ( _pid == 0 )
                std::cerr << "ABORTING PROGRAM:\n\tPlease, use a number of processes equal to NP*NP\n";
            abort( );
        }
        if ( _gsize%_dim != 0 )
        {
            if ( _pid == 0 )
                std::cerr << "ABORTING PROGRAM:\n\tSize of the raster [" << _gsize << "] is not multiple of N [" << _dim << "]\n";
            abort( );
        }
        _lsize = _gsize/_dim;
        if ( _lsize%2 != 0 )
        {
            if ( _pid == 0 )
                std::cerr << "ABORTING PROGRAM:\n\tLocal size of the raster [" << _lsize << "] is not multiple of 2\n";
            abort( );
        }

        init( );
    }

    OverlapAreas::~OverlapAreas( )
    {
    }

    void OverlapAreas::init( )
    {
        _worldPos = Point2D<int>( _pid%_dim, _pid/_dim );

        _area = Rectangle<int>( Size<int>( _lsize, _lsize ), Point2D<int>( _worldPos._x*_lsize, _worldPos._y*_lsize ) );
        _inn_area = Rectangle<int>( Size<int>( _lsize - 2*_overlap, _lsize - 2*_overlap ),
                                    Point2D<int>( _worldPos._x*_lsize  + _overlap, _worldPos._y*_lsize + _overlap ) );
        _ext_area = Rectangle<int>( std::max( _worldPos._x*_lsize - _overlap, 0 ),
                                    std::max( _worldPos._y*_lsize - _overlap, 0 ),
                                    std::min( (_worldPos._x+1)*_lsize - 1 + _overlap, _gsize - 1 ),
                                    std::min( (_worldPos._y+1)*_lsize - 1 + _overlap, _gsize - 1 ) );
        _sectionArea[0] = Rectangle<int>( Size<int>( _lsize/2, _lsize/2 ), Point2D<int>( _worldPos._x*_lsize, _worldPos._y*_lsize ) );
        _sectionArea[1] = _sectionArea[0] + Point2D<int>( _lsize/2, 0 );
        _sectionArea[2] = _sectionArea[0] + Point2D<int>( 0, _lsize/2 );
        _sectionArea[3] = _sectionArea[0] + Point2D<int>( _lsize/2, _lsize/2 );

        Interval<int> limits = Interval<int>( 0, _dim-1 );
        std::vector<Point2D<int>> increments =
        {
            Point2D<int>(-1,0), Point2D<int>(1,0), Point2D<int>(0,-1), Point2D<int>(0,1),
            Point2D<int>(-1,-1), Point2D<int>(1,1), Point2D<int>(1,-1), Point2D<int>(-1,1)
        };

        int ev = _worldPos._x%2 != _worldPos._y%2;
        for (int i= 0; i< 8; i++ )
        {
            int n = i+ev;
            Point2D<int> coord = _worldPos + increments[n];
            if ( limits.isInside( coord._x ) && limits.isInside( coord._y ) )
            {
                _neighbors.push_back(
                    std::make_tuple(
                        coord._x + _dim*coord._y,
                        Rectangle<int>( Size<int>( _lsize, _lsize ), Point2D<int>( coord._x*_lsize, coord._y*_lsize ) ),
                        Rectangle<int>( std::max( coord._x*_lsize - _overlap, 0 ),
                                        std::max( coord._y*_lsize - _overlap, 0 ),
                                        std::min( (coord._x+1)*_lsize - 1 + _overlap, _gsize - 1 ),
                                        std::min( (coord._y+1)*_lsize - 1 + _overlap, _gsize - 1 ) ) ) );
            }
            ev = -ev;
        }
    }

    const Rectangle<int> & OverlapAreas::getOwnedArea( ) const
    {
        return _area;
    }

    const Rectangle<int> & OverlapAreas::getInnerArea( ) const
    {
        return _inn_area;
    }

    const Rectangle<int> & OverlapAreas::getOverlapArea( ) const
    {
        return _ext_area;
    }

    Point2D<int> & OverlapAreas::getWorldPosition( )
    {
        return _worldPos;
    }

    Point2D<int> OverlapAreas::getRandomPosition( ) const
    {
        Engine::Point2D<int> pos( GeneralState::statistics( ).getUniformDistValue( 0, _lsize - 1 ),
                                  GeneralState::statistics( ).getUniformDistValue( 0, _lsize - 1 ) );
        return _area._origin + pos;
    }

    Point2D<int> OverlapAreas::getWorldOrigin( )
    {
        return _ext_area._origin;
    }

    NeightborVector &OverlapAreas::getNeighbors( )
    {
        return _neighbors;
    }

    int OverlapAreas::getSection( Point2D<int> pos ) const
    {
        if ( _sectionArea[0].contains( pos ) ) return 0;
        if ( _sectionArea[1].contains( pos ) ) return 1;
        if ( _sectionArea[2].contains( pos ) ) return 2;
        if ( _sectionArea[3].contains( pos ) ) return 3;
        return -1;
    }
#endif
}
