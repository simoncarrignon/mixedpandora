#ifndef __OverlapAreas_hxx__
#define __OverlapAreas_hxx__

#include <typedefs.hxx>
#include <Rectangle.hxx>
#include <Point2D.hxx>
#include <GeneralState.hxx>
#include <vector>
#include <tuple>

namespace Engine
{
    typedef std::vector<std::tuple<int, Rectangle<int>, Rectangle<int>>> NeightborVector;
    class OverlapAreas
    {
    private :
        int _pid;
        int _nprocs;
        int _overlap;
        int _gsize;
        int _dim;
        int _lsize;

        Point2D<int> _worldPos;

        Rectangle<int> _area;
        Rectangle<int> _ext_area;

        NeightborVector _neighbors;

        void init( );

    public :
        OverlapAreas( ) { };
        OverlapAreas( int pid, int nprocs, int overlap, int gsize );
        ~OverlapAreas( );

        const Rectangle<int> & getOwnedArea( ) const;
        const Rectangle<int> & getOverlapArea( ) const;
        Point2D<int> & getWorldPosition( );
        Point2D<int> getRandomPosition( ) const;
        Point2D<int> getWorldOrigin( );
        NeightborVector &getNeighbors( );

    }; // class OverlapAreas
}
#endif
