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

    typedef struct OverlapArea
    {
        int _n;
        Rectangle<int> _local;
        Rectangle<int> _bound;
    } Overlap_st;

    class OverlapAreas
    {
    private :
        int _pid;
        int _dim;

        Point2D<int> _pos;
        Rectangle<int> _area;
        Rectangle<int> _boun;

        Overlap_st _left;
        Overlap_st _topleft;
        Overlap_st _topright;

        int _l_neig;
        Rectangle<int> _l_area;
        Rectangle<int> _l_boun;

        Overlap_st _right;
        int _r_neig;
        Rectangle<int> _r_area;
        Rectangle<int> _r_boun;

        Overlap_st _top;
        int _t_neig;
        Rectangle<int> _t_area;
        Rectangle<int> _t_boun;

        Overlap_st _bottom;
        int _b_neig;
        Rectangle<int> _b_area;
        Rectangle<int> _b_boun;

        std::vector<int> _v_neig;
        std::vector<Rectangle<int>> _v_area;
        std::vector<Rectangle<int>> _v_boun;

        void push_back( int neig, Rectangle<int> area, Rectangle<int> boun );
    public :
        OverlapAreas( int id, int nprocs, int size, int over );
        void abort( );
        bool isEven( ) { return _pos._x%2 == _pos._y%2; }

        Overlap_st *getLeft( ) { return _pos._x > 0 ? &_left : 0; }
        Overlap_st *getRight( ) { return _pos._x < _dim-1 ? &_right : 0; }
        Overlap_st *getTop( ) { return _pos._y >0 ? &_top : 0; }
        Overlap_st *getBottom( ) { return _pos._y < _dim-1 ? &_bottom : 0; }

        Overlap_st *getTopLeft( ) { return _pos._x > 0 ? &_topleft : 0; }
        Overlap_st *getTopRight( ) { return _pos._x < _dim-1 ? &_topright : 0; }

        bool isLeft( ) { return _pos._x > 0; }
        bool isRight( ) { return _pos._x < _dim-1; }
        bool isTop( ) { return _pos._y > 0; }
        bool isBottom( ) { return _pos._y < _dim-1; }

        int getNumOfBounds( ) { return _v_neig.size(); }
        Rectangle<int> &getArea( ) { return _area; }
        Rectangle<int> &getBoundaries( ) { return _boun; }
        int getNeightbour( int i ) { return _v_neig[i]; }
        Rectangle<int> &getArea( int i ) { return _v_area[i]; }
        Rectangle<int> &getBoundary( int i ) { return _v_boun[i]; }
        std::vector<Rectangle<int>> &getAreas( ) { return _v_area; }
#ifdef RGT
        int _pid;
        int _nprocs;
        int _overlap;
        int _gsize;
        int _dim;
        int _lsize;

        Point2D<int> _worldPos;

        Rectangle<int> _area;
        Rectangle<int> _ext_area;
        Rectangle<int> _inn_area;

        Rectangle<int> _sectionArea[4];

        NeightborVector _neighbors;

        void init( );

    public :
        OverlapAreas( ) { };
        OverlapAreas( int pid, int nprocs, int overlap, int gsize );
        ~OverlapAreas( );

        const Rectangle<int> & getOwnedArea( ) const;
        const Rectangle<int> & getInnerArea( ) const;
        const Rectangle<int> & getOverlapArea( ) const;
        const Rectangle<int> & getSectionArea( int i ) const { return _sectionArea[i]; };
        Point2D<int> & getWorldPosition( );
        Point2D<int> getRandomPosition( ) const;
        Point2D<int> getWorldOrigin( );
        NeightborVector &getNeighbors( );

        int getSection( Point2D<int> pos ) const ;
#endif
    }; // class OverlapAreas
}
#endif
