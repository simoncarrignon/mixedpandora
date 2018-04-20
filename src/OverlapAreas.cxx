#include <OverlapAreas.hxx>
#include <Size.hxx>
#include <iostream>

namespace Engine
{
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
}
