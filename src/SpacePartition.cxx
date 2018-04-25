/*
 * Copyright ( c ) 2014
 * COMPUTER APPLICATIONS IN SCIENCE & ENGINEERING
 * BARCELONA SUPERCOMPUTING CENTRE - CENTRO NACIONAL DE SUPERCOMPUTACIÓN
 * http://www.bsc.es

 * This file is part of Pandora Library. This library is free software;
 * you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation;
 * either version 3.0 of the License, or ( at your option ) any later version.
 *
 * Pandora is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <SpacePartition.hxx>
#include <Agent.hxx>
#include <MpiFactory.hxx>
#include <Logger.hxx>
#include <Exception.hxx>
#include <Config.hxx>

namespace Engine
{
    SpacePartition::SpacePartition( const int & overlap, bool finalize ) : _serializer( *this ), _worldPos( -1, -1 ), _overlap( overlap ), _finalize( finalize ), _initialTime( 0.0f )
    {
    }

    SpacePartition::~SpacePartition( )
    {
        std::cout << "SpacePartition::~SpacePartition\n"; exit(1);
    }

    void SpacePartition::init( int argc, char *argv[] )
    {
        int alreadyInitialized;
        MPI_Initialized( &alreadyInitialized );
        if ( !alreadyInitialized )
        {
            MPI_Init( &argc, &argv );
        }
        _initialTime = getWallTime( );

        MPI_Comm_size( MPI_COMM_WORLD, &_numTasks );
        MPI_Comm_rank( MPI_COMM_WORLD, &_id );
        stablishBoundaries( );
    }

    void SpacePartition::initData( )
    {
        // serializer init
        _serializer.init( *_world );

        // mpi type registering
        MpiFactory::instance( )->registerTypes( );
        initOverlappingData( );

        std::stringstream logName;
        logName << "simulation_" << _id;
        log_INFO( logName.str( ), "finished init at: "  << getWallTime( ) );
    }

    void SpacePartition::checkOverlapSize( )
    {
        std::cout << "SpacePartition::checkOverlapSize\n"; exit(1);
#ifdef RGT
        int subfieldSizeX = _ownedArea._size._width/2;
        int subfieldSizeY = _ownedArea._size._height/2;
        if ( _overlap*2>subfieldSizeX || _overlap*2>subfieldSizeY )
        {
            std::stringstream oss;
            oss << "SpacePartition::checkOverlapSize- subfield sizes: " << subfieldSizeX << "/" << subfieldSizeY << " from global: " << _world->getConfig( ).getSize( ) << " and owned area: " << _ownedArea << " must be at least twice the value of overlap: " << _overlap << " to avoid conflicts between non adjacent subfields";
            throw Exception( oss.str( ) );
        }
#endif
    }

    void SpacePartition::stablishBoundaries( )
    {
        _mpiOverlap = new OverlapAreas( _id, _numTasks, _world->getConfig( ).getSize( )._width, _overlap );
        _ownedArea = _mpiOverlap->getArea( );
        _boundaries = _mpiOverlap->getBoundaries( );
#ifdef RGT
        int dim = sqrt( _numTasks );
        if ( dim*dim != _numTasks )
        {
            if ( _id == 0 )
                std::cerr << "ABORTING PROGRAM:\n\tPlease, use a number of processes equal to NP*NP\n";
            abort( );
        }
        int gsize = _world->getConfig( ).getSize( )._width;
        if ( gsize%dim != 0 )
        {
            if ( _id == 0 )
                std::cerr << "ABORTING PROGRAM:\n\tSize of the raster [" << gsize << "] is not multiple of N [" << dim << "]\n";
            abort( );
        }
        int lsize = gsize/dim;
        if ( lsize%2 != 0 )
        {
            if ( _id == 0 )
                std::cerr << "ABORTING PROGRAM:\n\tLocal size of the raster [" << lsize << "] is not multiple of 2\n";
            abort( );
        }
        Rectangle<int> global = Rectangle<int>( 0, 0, gsize-1, gsize-1 );
        _worldPos = Point2D<int>( _id%dim, _id/dim );

        _ownedArea = Rectangle<int>( Size<int>( lsize, lsize ), Point2D<int>( _worldPos._x*lsize, _worldPos._y*lsize ) );
        _boundaries = global.intersection( _ownedArea - Point2D<int>( _overlap, _overlap ) + Size<int>( 2*_overlap, 2*_overlap ) );

        // LEFT BORDER
        int ind = _worldPos._x%2 == _worldPos._y%2 ? 0 : 1;
        if ( _worldPos._x > 0 )
        {
            _nearby[ind] = _id - 1;
            _nb_area[ind] = Rectangle<int>( _ownedArea.left( ), _ownedArea.top( ), _ownedArea.left( )+_overlap-1, _ownedArea.bottom() );
            _nb_exta[ind] = Rectangle<int>( _ownedArea.left( )-_overlap, _ownedArea.top( ), _ownedArea.left( )-1, _ownedArea.bottom() );
        }
        else
            _nearby[ind] = -1;

        // RIGHT BORDER
        ind = _worldPos._x%2 == _worldPos._y%2 ? 1 : 0;
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

        _left_top = _worldPos._x > 0 ? Rectangle<int>( _boundaries.left( ), _boundaries.top( ),
                                       _ownedArea.left( )+_overlap-1, _ownedArea.top( )+lsize/2 + _overlap - 1 ) : Rectangle<int>( );
        _left_bot = _worldPos._x > 0 ? Rectangle<int>( _boundaries.left( ), _ownedArea.top( ) + lsize/2 - _overlap,
                                        _ownedArea.left( )+ _overlap - 1, _boundaries.bottom() ) :  Rectangle<int>( );
        _top_broad = _worldPos._y > 0 ? Rectangle<int>( _boundaries.left( ), _boundaries.top( ),
                                                        _boundaries.right( ), _ownedArea.top( ) + _overlap - 1 ) : Rectangle<int>( );
        _top_bou = _worldPos._y > 0 ? Rectangle<int>( _ownedArea.left( ), _boundaries.top( ),
                                      _ownedArea.right( ), _ownedArea.top( ) + _overlap - 1 ) : Rectangle<int>( );
        _right_top = _worldPos._x < dim-1 ? Rectangle<int>( _ownedArea.right( ) - _overlap + 1, _boundaries.top( ),
                                            _boundaries.right( ), _ownedArea.top( ) + lsize/2 + _overlap - 1 ) : Rectangle<int>( );
//        if ( _neighbors[0] )
//            _left_bound = 

#endif

//        _inn_area = _OverlapAreas.getInnerArea( );
//        _boundaries = _OverlapAreas.getOverlapArea( );
#ifdef KKK
        // position of world related to the complete set of computer nodes
        int worldsPerRow = sqrt( _numTasks );
        int size = _world->getConfig( ).getSize( )._width;
        size = size/worldsPerRow;
        _worldPos = getPositionFromId( _id );
#ifdef ORIG
        for ( int x=_worldPos._x-1; x<=_worldPos._x+1; x++ )
        {
            for ( int y=_worldPos._y-1; y<=_worldPos._y+1; y++ )
            {
                if ( x>-1 && x<worldsPerRow && y>-1 && y<worldsPerRow )
                {
                    if ( x!=_worldPos._x || y!=_worldPos._y )
                    {
                        _neighbors.push_back( y*worldsPerRow+x );
                    }
                }
            }
        }
        // owned area inside global coordinates, depending on worldPos
        _ownedArea._size._width = _world->getConfig( ).getSize( )._width/worldsPerRow;
        _ownedArea._size._height = _world->getConfig( ).getSize( )._height/worldsPerRow;
        _ownedArea._origin._x = _worldPos._x*_ownedArea._size._width;
        _ownedArea._origin._y = _worldPos._y*_ownedArea._size._height;

        // defining overlap boundaries
        _boundaries = _ownedArea;
        // west boundary
        if ( _ownedArea._origin._x!=0 )
        {
            _boundaries._origin._x -= _overlap;
            _boundaries._size._width += _overlap;
        }
        // east boundary
        if ( _ownedArea._origin._x!=_world->getConfig( ).getSize( )._width-_ownedArea._size._width )
        {
            _boundaries._size._width += _overlap;
        }
        // north boundary
        if ( _ownedArea._origin._y!=0 )
        {
            _boundaries._origin._y -= _overlap;
            _boundaries._size._height += _overlap;
        }
        // south boundary
        if ( _ownedArea._origin._y!=_world->getConfig( ).getSize( )._height-_ownedArea._size._height )
        {
            _boundaries._size._height += _overlap;
        }

        if ( _ownedArea._size._width%2!=0 || _ownedArea._size._height%2!=0 )
        {
            std::stringstream oss;
            oss << "SpacePartition::init - local raster size: " << _ownedArea._size << " must be divisible by 2";
            throw Exception( oss.str( ) );
        }
        std::cout << _id << " _ownedArea: " << _ownedArea << " " << _OverlapAreas.getOwnedArea( ) << "\n";
        std::cout << _id << " _boundaries: " << _boundaries << " " << _OverlapAreas.getOverlapArea( ) << "\n"; exit(0);

        checkOverlapSize( );

        // creating sections
        _sections.resize( 4 );
        _sections[0] = Rectangle<int>( _ownedArea._size/2, _ownedArea._origin );
        _sections[1] = _sections[0] + Point2D<int>( size/2, 0 );
        _sections[2] = _sections[0] + Point2D<int>( 0, size/2 );
        _sections[3] = _sections[0] + Point2D<int>( size/2, size/2 );

#else
#ifdef RGT
        setSectionNeighbours( worldsPerRow );
        // owned area inside global coordinates, depending on worldPos
        _ownedArea = Rectangle<int>( _worldPos._x*size, _worldPos._y*size,
                              ( _worldPos._x + 1 )*size-1, ( _worldPos._y + 1 )*size-1 );
        
        _boundaries = Rectangle<int>( std::max( 0, _worldPos._x*size - _overlap ), std::max( 0, _worldPos._y*size - _overlap ),
            std::min( size*worldsPerRow, ( _worldPos._x + 1 )*size + _overlap ) - 1, std::min( size*worldsPerRow, ( _worldPos._y + 1 )*size + _overlap ) - 1 );

        // creating sections
        _sections.resize( 4 );
        _sections[0] = Rectangle<int>( _ownedArea._size/2, _ownedArea._origin );
        _sections[1] = _sections[0] + Point2D<int>( size/2, 0 );
        _sections[2] = _sections[0] + Point2D<int>( 0, size/2 );
        _sections[3] = _sections[0] + Point2D<int>( size/2, size/2 );
#endif
#endif

        std::stringstream logName;
        logName << "simulation_" << _id;
        log_INFO( logName.str( ), getWallTime( ) << " pos: " << _worldPos << ", global size: " << _world->getConfig( ).getSize( ) << ", boundaries: " << _boundaries << " and owned area: " << _ownedArea );
        log_INFO( logName.str( ), getWallTime( ) << " sections 0: " << _sections[0] << " - 1: " << _sections[1] << " - 2:" << _sections[2] << " - 3: " << _sections[3] );
#endif
    }

#ifdef RGT
    void SpacePartition::setSectionNeighbours( int nc )
    {
        if (_worldPos._x%2==_worldPos._y)
        {
            // section 0
            if ( _worldPos._x>0 ) _Neighbors[0].push_back( std::make_pair( true, _id-1 ) );
            if ( _worldPos._y>0 ) _Neighbors[0].push_back( std::make_pair( true, _id-nc ) );
            if ( _worldPos._x>0 && _worldPos._y>0 ) _Neighbors[0].push_back( std::make_pair( true, _id-nc-1 ) );
            if ( _worldPos._x<nc-1 ) _Neighbors[0].push_back( std::make_pair( false, _id+1 ) );
            if ( _worldPos._y<nc-1 ) _Neighbors[0].push_back( std::make_pair( false, _id+nc ) );
            if ( _worldPos._x<nc-1 && _worldPos._y<nc-1 ) _Neighbors[0].push_back( std::make_pair( false, _id+1+nc ) );
            // section 1
            if ( _worldPos._x<nc-1 ) _Neighbors[1].push_back( std::make_pair( true, _id+1 ) );
            if ( _worldPos._y>0 ) _Neighbors[1].push_back( std::make_pair( true, _id-nc ) );
            if ( _worldPos._x<nc-1 && _worldPos._y>0 ) _Neighbors[1].push_back( std::make_pair( true, _id-nc+1 ) );
            if ( _worldPos._x>0 ) _Neighbors[1].push_back( std::make_pair( false, _id-1 ) );
            if ( _worldPos._y<nc-1 ) _Neighbors[1].push_back( std::make_pair( false, _id+nc ) );
            if ( _worldPos._x>0 && _worldPos._y<nc-1 ) _Neighbors[1].push_back( std::make_pair( false, _id-1+nc ) );
            // section 2
            if ( _worldPos._x>0 ) _Neighbors[2].push_back( std::make_pair( true, _id-1 ) );
            if ( _worldPos._y<nc-1 ) _Neighbors[2].push_back( std::make_pair( true, _id+nc ) );
            if ( _worldPos._x>0 && _worldPos._y<nc-1 ) _Neighbors[2].push_back( std::make_pair( true, _id+nc-1 ) );
            if ( _worldPos._x<nc-1 ) _Neighbors[2].push_back( std::make_pair( false, _id+1 ) );
            if ( _worldPos._y>0 ) _Neighbors[2].push_back( std::make_pair( false, _id-nc ) );
            if ( _worldPos._x<nc-1 && _worldPos._y<nc-1 ) _Neighbors[2].push_back( std::make_pair( false, _id+1-nc ) );
            // section 3
            if ( _worldPos._x<nc-1 ) _Neighbors[1].push_back( std::make_pair( true, _id+1 ) );
            if ( _worldPos._y<nc-1 ) _Neighbors[1].push_back( std::make_pair( true, _id+nc ) );
            if ( _worldPos._x<nc-1 && _worldPos._y>nc-1 ) _Neighbors[1].push_back( std::make_pair( true, _id+nc+1 ) );
            if ( _worldPos._x>0 ) _Neighbors[1].push_back( std::make_pair( false, _id-1 ) );
            if ( _worldPos._y>0 ) _Neighbors[1].push_back( std::make_pair( false, _id-nc ) );
            if ( _worldPos._x>0 && _worldPos._y>0 ) _Neighbors[1].push_back( std::make_pair( false, _id-1-nc ) );
        }
        else
        {
            // section 0
            if ( _worldPos._x<nc-1 ) _Neighbors[0].push_back( std::make_pair( false, _id+1 ) );
            if ( _worldPos._y<nc-1 ) _Neighbors[0].push_back( std::make_pair( false, _id+nc ) );
            if ( _worldPos._x<nc-1 && _worldPos._y<nc-1 ) _Neighbors[0].push_back( std::make_pair( false, _id+1+nc ) );
            if ( _worldPos._x>0 ) _Neighbors[0].push_back( std::make_pair( true, _id-1 ) );
            if ( _worldPos._y>0 ) _Neighbors[0].push_back( std::make_pair( true, _id-nc ) );
            if ( _worldPos._x>0 && _worldPos._y>0 ) _Neighbors[0].push_back( std::make_pair( true, _id-nc-1 ) );
            // section 1
            if ( _worldPos._x>0 ) _Neighbors[1].push_back( std::make_pair( false, _id-1 ) );
            if ( _worldPos._y<nc-1 ) _Neighbors[1].push_back( std::make_pair( false, _id+nc ) );
            if ( _worldPos._x>0  && _worldPos._y<nc-1 ) _Neighbors[1].push_back( std::make_pair( false, _id-1+nc ) );
            if ( _worldPos._x<nc-1 ) _Neighbors[1].push_back( std::make_pair( true, _id+1 ) );
            if ( _worldPos._y>0 ) _Neighbors[1].push_back( std::make_pair( true, _id-nc ) );
            if ( _worldPos._x<nc-1 && _worldPos._y>0 ) _Neighbors[1].push_back( std::make_pair( true, _id-nc+1 ) );
            // section 2
            if ( _worldPos._x<nc-1 ) _Neighbors[2].push_back( std::make_pair( false, _id+1 ) );
            if ( _worldPos._y>0 ) _Neighbors[2].push_back( std::make_pair( false, _id-nc ) );
            if ( _worldPos._x<nc-1 && _worldPos._y<nc-1 ) _Neighbors[2].push_back( std::make_pair( false, _id+1-nc ) );
            if ( _worldPos._x>0 ) _Neighbors[2].push_back( std::make_pair( true, _id-1 ) );
            if ( _worldPos._y<nc-1 ) _Neighbors[2].push_back( std::make_pair( true, _id+nc ) );
            if ( _worldPos._x>0 && _worldPos._y<nc-1 ) _Neighbors[2].push_back( std::make_pair( true, _id+nc-1 ) );
            // section 3
            if ( _worldPos._x>0 ) _Neighbors[3].push_back( std::make_pair( false, _id-1 ) );
            if ( _worldPos._y>0 ) _Neighbors[3].push_back( std::make_pair( false, _id-nc ) );
            if ( _worldPos._x>0 && _worldPos._y>0 ) _Neighbors[3].push_back( std::make_pair( false, _id-1-nc ) );
            if ( _worldPos._x<nc-1 ) _Neighbors[3].push_back( std::make_pair( true, _id+1 ) );
            if ( _worldPos._y<nc-1 ) _Neighbors[3].push_back( std::make_pair( true, _id+nc ) );
            if ( _worldPos._x<nc-1 && _worldPos._y>nc-1 ) _Neighbors[3].push_back( std::make_pair( true, _id+nc+1 ) );

        }
    }
#else
#endif

#ifdef ORIG
    void SpacePartition::stepSection( const int & sectionIndex )
    {
        std::cout << "SpacePartition::stepSection\n"; MPI_Barrier( MPI_COMM_WORLD ); MPI_Finalize( ); exit(0);
        std::stringstream logName;
        logName << "simulation_" << _id;
        log_DEBUG( logName.str( ), getWallTime( ) << " beginning step: " << _world->getCurrentStep( ) << " section: " << sectionIndex );

        AgentsList::iterator it = _world->beginAgents( );
        AgentsVector agentsToExecute;
        // we have to randomize the execution of agents in a given section index
        while( it!=_world->endAgents( ) )
        {
            AgentPtr agent = *it;
            if ( _sections[sectionIndex].contains( agent->getPosition( ) ) && !hasBeenExecuted( agent->getId( ) ))
            {
                log_DEBUG( logName.str( ), "agent to execute: :" << agent );
                agentsToExecute.push_back( agent );
            }
            it++;
        }
        std::random_shuffle( agentsToExecute.begin( ), agentsToExecute.end( ) );
        int numExecutedAgents = 0;
        AgentsList agentsToSend;


    #ifndef PANDORAEDEBUG
        // shared memory distibution for read-only planning actions, disabled for extreme debug
        #pragma omp parallel for
    #endif
        for ( size_t i=0; i<agentsToExecute.size( ); i++ )
        {
            agentsToExecute.at( i )->updateKnowledge( );
            agentsToExecute.at( i )->selectActions( );
        }

        // execute actions
        for ( size_t i=0; i<agentsToExecute.size( ); i++ )
        {
            AgentPtr agent = agentsToExecute.at( i );
            log_DEBUG( logName.str( ), getWallTime( ) << " agent: " << agent << " being executed at index: " << sectionIndex << " of task: "<< _id << " in step: " << _world->getCurrentStep( ) );

            agentsToExecute.at( i )->executeActions( );
            agentsToExecute.at( i )->updateState( );
            log_DEBUG( logName.str( ), getWallTime( ) << " agent: " << agent << " has been executed at index: " << sectionIndex << " of task: "<< _id << " in step: " << _world->getCurrentStep( ) );

            if ( !_ownedArea.contains( agent->getPosition( ) ) && !willBeRemoved( agent->getId( ) ))
            {
                log_DEBUG( logName.str( ), getWallTime( ) << " migrating agent: " << agent << " being executed at index: " << sectionIndex << " of task: "<< _id );
                std::cout << _id << "  Migrating agent: " << agent << "\n";
                agentsToSend.push_back( agent );

                // the agent is no longer property of this world
                AgentsList::iterator itErase  = getOwnedAgentPtr( agent->getId( ) );
                // it will be deleted
                _world->eraseAgent( itErase );
                _overlapAgents.push_back( agent );

                std::cout << _id << "  ag_T_Send: " << agentsToSend.size() << "  n_ag: " << _world->getNumberOfAgents()  << "  Ov_ag: " << _overlapAgents.size()<< "\n";

                log_DEBUG( logName.str( ), getWallTime( ) <<  "putting agent: " << agent << " to overlap" );
            }
            else
            {
                log_DEBUG( logName.str( ), getWallTime( ) << " finished agent: " << agent );
            }
            _executedAgentsHash.insert( make_pair( agent->getId( ), agent ));
            numExecutedAgents++;
            log_DEBUG( logName.str( ), getWallTime( )  << " num executed agents: " << numExecutedAgents );
        }
        log_DEBUG( logName.str( ), getWallTime( )  << " sending agents in section: " << sectionIndex << " and step: " << _world->getCurrentStep( ) );
        sendAgents( agentsToSend );
        log_DEBUG( logName.str( ), getWallTime( ) << " has finished section: " << sectionIndex << " and step: " << _world->getCurrentStep( ) );

        log_DEBUG( logName.str( ), getWallTime( ) << " executed step: " << _world->getCurrentStep( ) << " section: " << sectionIndex << " in zone: " << _sections[sectionIndex] << " with num executed agents: " << numExecutedAgents << " total agents: " << std::distance( _world->beginAgents( ), _world->endAgents( ) ) << " and overlap agents: " << _overlapAgents.size( ) );
    }
#else
    void SpacePartition::stepSection( const int & sectionIndex, AgentsVector &agentsToExecute )
    {
#ifdef TMP
        std::random_shuffle( agentsToExecute.begin( ), agentsToExecute.end( ) );
        std::cout << _id << " " << sectionIndex << " SpacePartition::stepSection " << agentsToExecute.size() << "\n";
        int numExecutedAgents = 0;
        AgentsList agentsToSend;
        for ( auto agent: agentsToExecute )
        {
            agent->updateKnowledge( );
            agent->selectActions( );
        }
        for ( auto agent: agentsToExecute )
        {
            agent->executeActions( );
            agent->updateState( );

            if ( !_ownedArea.contains( agent->getPosition( ) ) )
            {
                std::cout << _id << " Migrating agent: " << agent->getId( ) << " to be removed? " << willBeRemoved( agent->getId( ) ) << "\n";
                agentsToSend.push_back( agent );
                migrateAgent( agent );
            }
            else if ( !_inn_area.contains( agent->getPosition( ) ) )
            {
                agentsToSend.push_back( agent );
                std::cout << _id << " Inner ghost agent " << agent->getId( ) << " to be removed? " << willBeRemoved( agent->getId( ) ) << "\n";
            }
        }
#endif   
    }
#endif

#ifdef TMP
    void SpacePartition::migrateAgent( AgentPtr agent )
    {
        for ( AgentsList::iterator it=_world->beginAgents( ); it!=_world->endAgents( ); it++ )
        {
            if ( (*it )->getId( ).compare( agent->getId( ) )==0 )
            {
                _world->eraseAgent( it );
                break;
            }
        }
        _overlapAgents.push_back( agent );
    }
#endif

    void SpacePartition::sendAgents( AgentsList & agentsToSend )
    {
        std::cout << "SpacePartition::sendAgents\n"; exit(1);
#ifdef ORIG
        if ( _neighbors.size( )==0 )
        {
            return;
        }
        std::stringstream logName;
        logName << "MPI_agents_world_" << _id;
        log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " sendAgent: " << agentsToSend.size( ) << " agents" );

        // for each neighbor, we send the number of agents to send
        for ( MpiFactory::TypesMap::iterator itType=MpiFactory::instance( )->beginTypes( ); itType!=MpiFactory::instance( )->endTypes( ); itType++ )
        {
            log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " sendAgent - checking mpi type: " << itType->first );
            // add each agent to the list of the neighbour where it will be sent
            std::vector< AgentsList > agentsToNeighbors;
            agentsToNeighbors.resize( _neighbors.size( ) );

            for ( AgentsList::iterator it=agentsToSend.begin( ); it!=agentsToSend.end( ); it++ )
            {
                AgentPtr agent = *it;
                if ( agent->isType( itType->first ))
                {
                    int newID = getIdFromPosition( agent->getPosition( ) );
                    agentsToNeighbors[getNeighborIndex( newID )].push_back( agent );
                }
            }

            MPI_Datatype * agentType = itType->second;
            for ( size_t i=0; i<_neighbors.size( ); i++ )
            {
                int numAgents = agentsToNeighbors[i].size( );
                log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " sendAgent - sending num agents: " << numAgents << " to: " << _neighbors[i] );

                int error = MPI_Send( &numAgents, 1, MPI_INTEGER, _neighbors[i], eNumAgents, MPI_COMM_WORLD );
                if ( error != MPI_SUCCESS )
                {
                    std::stringstream oss;
                    oss << "SpacePartition::sendAgents - " << _id << " error in MPI_Send : " << error;
                    throw Exception( oss.str( ) );
                }
                AgentsList::iterator it= agentsToNeighbors[i].begin( );
                while( it!=agentsToNeighbors[i].end( ) )
                {
                    Agent * agent = it->get( );
                    void * package = agent->fillPackage( );
                    log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " sendAgent - sending agent: " << *it << " to: " << _neighbors[i] );
                    error = MPI_Send( package, 1, *agentType, _neighbors[i], eAgent, MPI_COMM_WORLD );
                    delete package;
                    if ( error != MPI_SUCCESS )
                    {
                        std::stringstream oss;
                        oss << "SpacePartition::sendAgents - " << _id << " error in MPI_Send : " << error;
                        throw Exception( oss.str( ) );
                    }
                    agent->sendVectorAttributes( _neighbors[i] );
                    // it is not deleted, as it is sent to overlap list
                    it = agentsToNeighbors[i].erase( it );
                }
            }
        }
        log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " sendAgent -  end checking agents to send: " << agentsToSend.size( ) );
#endif
    }

    void SpacePartition::sendOverlapZones( const int & sectionIndex, const bool & entireOverlap )
    {
        std::cout << "SpacePartition::sendOverlapZones\n"; exit(1);
#ifdef ORIG
        std::stringstream logName;
        logName << "MPI_raster_world_" << _id;
        log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << "/" << sectionIndex << " sendOverlapZones" );
        std::vector<int> neighborsToUpdate;

        for ( size_t i=0; i<_neighbors.size( ); i++ )
        {
            if ( needsToBeUpdated( _neighbors[i], sectionIndex ))
            {
                neighborsToUpdate.push_back( _neighbors[i] );
            }
        }

        for ( size_t d=0; d<_world->getNumberOfRasters( ); d++ )
        {
            if ( !_world->rasterExists( d ) || !_world->isRasterDynamic( d ))
            {
                log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " index: " << d << " not sending" );
                continue;
            }
            log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << "/" << sectionIndex << " sending raster: " << d );
            for ( size_t i=0; i<neighborsToUpdate.size( ); i++ )
            {
                MpiOverlap * send = new MpiOverlap;
                if ( entireOverlap )
                {
                    send->_overlap= getOverlap( neighborsToUpdate[i], sectionIndex );
                    log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << "/" << sectionIndex << " send entire overlap: " << send->_overlap << " to: " << neighborsToUpdate[i] );
                }
                else
                {
                    send->_overlap = getInternalOverlap( neighborsToUpdate[i] );
                    log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << "/" << sectionIndex << " send partial overlap: " << send->_overlap << " to: " << neighborsToUpdate[i] );
                }
                const Rectangle<int> & overlapZone = send->_overlap;
                send->_data.resize( overlapZone._size._width * overlapZone._size._height );
                log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << "/" << sectionIndex << " will send overlap to: " << neighborsToUpdate[i] << " with size: " << send->_data.size( ) << " and zone: " << overlapZone );
                for ( size_t n=0; n<send->_data.size( ); n++ )
                {
                    Point2D<int> index( overlapZone._origin._x+n%overlapZone._size._width, overlapZone._origin._y+n/overlapZone._size._width );
                    send->_data.at( n ) = _world->getDynamicRaster( d ).getValue( index );
                    log_EDEBUG( logName.str( ), "\t" << getWallTime( ) << " step: " << _world->getCurrentStep( ) << "/" << sectionIndex << " send index: " << index << " in global pos: " << index+_boundaries._origin << " value: " << send->_data.at( n ));
                }
                log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " raster: " << d << " will be sent" );
                MPI_Isend( &send->_data[0], send->_data.size( ), MPI_INTEGER, neighborsToUpdate[i], eRasterData, MPI_COMM_WORLD, &send->_request );
                _sendRequests.push_back( send );
                log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " raster: " << d << " data sent to: " << _neighbors[i] );
            }
            log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << "/" << sectionIndex << " raster: " << d << " sent" );
        }
        log_DEBUG( logName.str( ), getWallTime( ) << " step: "  << "/" << sectionIndex << _world->getCurrentStep( ) << " sendOverlapZones ended" );
#endif
    }

    void SpacePartition::sendMaxOverlapZones( )
    {
#ifdef ORIG
        std::stringstream logName;
        logName << "MPI_raster_world_" << _id;
        log_DEBUG( logName.str( ), getWallTime( ) << " step: "  << _world->getCurrentStep( ) << " sendMaxOverlapZones" );
        for ( size_t d=0; d<_world->getNumberOfRasters( ); d++ )
        {
            if ( !_world->rasterExists( d ) || !_world->isRasterDynamic( d ))
            {
                continue;
            }

            log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " sending max raster: " << d );
            for ( size_t i=0; i<_neighbors.size( ); i++ )
            {
                MpiOverlap * send = new MpiOverlap;
                send->_overlap = getInternalOverlap( _neighbors[i] );
                send->_data.resize( send->_overlap._size._width * send->_overlap._size._height );
                log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " will send max overlap of: " << d << " to: " << _neighbors[i] << " with size: " << send->_data.size( ) << " and zone: " << send->_overlap << " to " << _neighbors[i] );
                const Rectangle<int> & overlapZone = send->_overlap;
                for ( size_t n=0; n<send->_data.size( ); n++ )
                {
                    Point2D<int> index( overlapZone._origin._x+n%overlapZone._size._width, overlapZone._origin._y+n/overlapZone._size._width );
                    send->_data.at( n ) = _world->getDynamicRaster( d ).getMaxValue( index );
                    log_EDEBUG( logName.str( ), "\t" << getWallTime( ) << " step: " << _world->getCurrentStep( ) << " send index: " << index << " in global pos: " << index+_boundaries._origin << " max value: " << send->_data.at( n ));
                }
                MPI_Isend( &send->_data[0], send->_data.size( ), MPI_INTEGER, _neighbors[i], eRasterMaxData, MPI_COMM_WORLD, &send->_request );
                _sendRequests.push_back( send );
                log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " raster: " << d << " max data sent to: " << _neighbors[i] );
            }
            log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " raster: " << d << " max data sent" );
        }
        log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " sendMaxOverlapZones ended" );
        std::cout << _id << " SpacePartition::sendMaxOverlapZones " << _world->getNumberOfRasters( ) << "\n"; exit(1);
#endif
    }

    void SpacePartition::sendOverlapZone( Rectangle<int> area, int dst )
    {
        std::vector<int> sdata( area.size()*2*_rasters.size() );
        area = area - _boundaries._origin;
        int i = 0;
        for ( auto raster: _rasters )
        {
            for ( auto coord: area )
            {
                sdata.at(i) = raster->getMaxValue( coord );
                sdata.at(i+1) = raster->getValue( coord );
                i += 2;
            }
        }
        MPI_Status status;
        MPI_Send( sdata.data( ), sdata.size(), MPI_INTEGER, dst, eRasterMaxData, MPI_COMM_WORLD );
    }

    void SpacePartition::receiveOverlapZone( Rectangle<int> area, int src )
    {
        std::vector<int> rdata( area.size()*2*_rasters.size() );
        MPI_Status status;
        MPI_Recv( rdata.data( ), rdata.size(), MPI_INTEGER, src, eRasterMaxData, MPI_COMM_WORLD, &status );
        area = area - _boundaries._origin;
        int i = 0;
        for ( auto raster: _rasters )
        {
            for ( auto coord: area )
            {
                raster->setMaxValue( coord, rdata[i] );
                raster->setValue( coord, rdata[i+1] );
                i+=2;
            }
        }
    }

    void SpacePartition::reduceOverlapZones( )
    {
        Overlap_st *b;
        // Send Right / Recv left
        if ( _mpiOverlap->isEven( ) )
        {
            if ( ( b = _mpiOverlap->getLeft( ) ) != 0 ) receiveOverlapZone( b->_bound, b->_n );
            if ( ( b = _mpiOverlap->getRight( ) ) != 0 ) sendOverlapZone( b->_local, b->_n );
            if ( ( b = _mpiOverlap->getTop( ) ) != 0 ) receiveOverlapZone( b->_bound, b->_n );
            if ( ( b = _mpiOverlap->getBottom( ) ) != 0 ) sendOverlapZone( b->_local, b->_n );
        }
        else
        {
            if ( ( b = _mpiOverlap->getRight( ) ) != 0 ) sendOverlapZone( b->_local, b->_n );
            if ( ( b = _mpiOverlap->getLeft( ) ) != 0 ) receiveOverlapZone( b->_bound, b->_n );
            if ( ( b = _mpiOverlap->getBottom( ) ) != 0 ) sendOverlapZone( b->_local, b->_n );
            if ( ( b = _mpiOverlap->getTop( ) ) != 0 ) sendOverlapZone( b->_local, b->_n );
        }
    }

    void SpacePartition::sendGhostAgents( const int & sectionIndex )
    {
        std::cout << "SpacePartition::sendGhostAgents\n"; exit(1);
#ifdef ORIG
        std::stringstream logName;
        logName << "MPI_agents_world_" << _id;
        log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " send ghost agents for section index: " << sectionIndex );

        std::vector<int> neighborsToUpdate;
        for ( size_t i=0; i<_neighbors.size( ); i++ )
        {
            if ( needsToBeUpdated( _neighbors[i], sectionIndex ))
            {
                neighborsToUpdate.push_back( _neighbors[i] );
                log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " section index: " << sectionIndex << " will send overlap to: " << _neighbors[i] );
            }
        }

        // for each type of agent we will send the collection of agents of the particular type to neighbors
        for ( MpiFactory::TypesMap::iterator itType=MpiFactory::instance( )->beginTypes( ); itType!=MpiFactory::instance( )->endTypes( ); itType++ )
        {
            log_DEBUG( logName.str( ),  getWallTime( ) << " step: " << _world->getCurrentStep( ) << " section index: " << sectionIndex << " checking type: " << itType->first );
            MPI_Datatype * agentType = itType->second;

            std::vector< AgentsList > agentsToNeighbors;
            agentsToNeighbors.resize( neighborsToUpdate.size( ) );

            for ( size_t i=0; i<neighborsToUpdate.size( ); i++ )
            {
                Rectangle<int> overlapZone = getOverlap( neighborsToUpdate[i], sectionIndex );
                for ( AgentsList::iterator it=_world->beginAgents( ); it!=_world->endAgents( ); it++ )
                {
                    AgentPtr agent = *it;
                    // we check the type. TODO register the type in another string
                    // TODO refactor!!!
                    log_DEBUG( logName.str( ),  getWallTime( ) << " step: " << _world->getCurrentStep( ) << " agent: " << agent << " of type: " << itType->first << " test will be removed: " << willBeRemoved( agent->getId( ) ) << " checking overlap zone: " << overlapZone << " overlap boundaries: " << _boundaries << " - test is inside zone: " << overlapZone.contains( agent->getPosition( )-_boundaries._origin ));
                    if ( agent->isType( itType->first ))
                    {
                        if ( (!willBeRemoved( agent->getId( ) )) && ( overlapZone.contains( agent->getPosition( )-_boundaries._origin )) )
                        {
                            agentsToNeighbors[i].push_back( agent );
                            log_DEBUG( logName.str( ),  getWallTime( ) << " step: " << _world->getCurrentStep( ) << " sending ghost agent: " << agent << " to: " << neighborsToUpdate[i] << " in section index: " << sectionIndex );
                        }
                    }
                }
                for ( AgentsList::iterator it=_overlapAgents.begin( ); it!=_overlapAgents.end( ); it++ )
                {
                    AgentPtr agent = *it;
                    if ( agent->isType( itType->first ))
                    {
                        if ( (!willBeRemoved( agent->getId( ) )) && ( overlapZone.contains( agent->getPosition( )-_boundaries._origin )) )
                        {
                            agentsToNeighbors[i].push_back( agent );
                            log_DEBUG( logName.str( ),  getWallTime( ) << " step: " << _world->getCurrentStep( ) << " will send modified ghost agent: " << agent << " to: " << neighborsToUpdate[i] << " in section index: " << sectionIndex << " and step: " << _world->getCurrentStep( ) );
                        }
                    }
                }

                int numAgents = agentsToNeighbors[i].size( );
                log_DEBUG( logName.str( ),  getWallTime( ) << " step: " << _world->getCurrentStep( ) << " sending num ghost agents: " << numAgents << " to : " << neighborsToUpdate[i] << " in step: " << _world->getCurrentStep( ) << " and section index: " << sectionIndex );
                int error = MPI_Send( &numAgents, 1, MPI_INTEGER, neighborsToUpdate[i], eNumGhostAgents, MPI_COMM_WORLD );
                if ( error != MPI_SUCCESS )
                {
                    std::stringstream oss;
                    oss << "SpacePartition::sendGhostAgents - " << _id << " error in MPI_Send : " << error;
                    throw Exception( oss.str( ) );
                }
                for ( AgentsList::iterator it=agentsToNeighbors[i].begin( ); it!=agentsToNeighbors[i].end( ); it++ )
                {
                    Agent * agent = it->get( );
                    void * package = agent->fillPackage( );
                    log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " sending ghost agent: " << *it << " from: " << _id << " to: " << neighborsToUpdate[i] );
                    error = MPI_Send( package, 1, *agentType, neighborsToUpdate[i], eGhostAgent, MPI_COMM_WORLD );
                    delete package;
                    if ( error != MPI_SUCCESS )
                    {
                        std::stringstream oss;
                        oss << "SpacePartition::sendGhostAgents - " << _id << " error in MPI_Send : " << error;
                        throw Exception( oss.str( ) );
                    }
                    agent->sendVectorAttributes( neighborsToUpdate[i] );
                }
            }
        }
        log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " send ghost agents for section index: " << sectionIndex << " finished" );
#endif
    }

    void SpacePartition::sendAgentsList( AgentsList list, int mpi_id, MPI_Datatype * agentType )
    {
        for ( AgentsList::iterator it=list.begin( ); it!=list.end( ); it++ )
        {
            Agent * agent = it->get( );
            void * package = agent->fillPackage( );
            int error = MPI_Send( package, 1, *agentType, mpi_id, eGhostAgent, MPI_COMM_WORLD );
            delete package;
            agent->sendVectorAttributes( mpi_id );
        }

    }

    void SpacePartition::receiveAgentsList( int nAgents, int mpi_id, const std::string & type, MPI_Datatype * agentType )
    {
        std::cout << _id << " receiveAgentsList\n"; exit(0);
#ifdef NEW
        for (int i=0; i< nAgents; i++ )
        {
            MPI_Status status;
            void * package = MpiFactory::instance( )->createDefaultPackage( type );
            int error = MPI_Recv( package, 1, *agentType, mpi_id, eGhostAgent, MPI_COMM_WORLD, &status );
            Agent * agent = MpiFactory::instance( )->createAndFillAgent( type, package );
            delete package;
            agent->receiveVectorAttributes( mpi_id );
            std::cout << _id << " " << agent->getId( ) << " -- " << mpi_id << " WORLD " << "\n";
            AgentsList::iterator it = getGhostAgentPtr( agent->getId( ) );
            if ( it != _overlapAgents.end( ) )
                _overlapAgents.erase( it );
            else
            {
                it = getOwnedAgentPtr( agent->getId( ) );
                if ( it!=_world->endAgents( ) )
                    _world->eraseAgent( it );
            }
            if ( _ownedArea.contains( agent->getPosition( ) ) )
                _world->addAgent( agent, false );
            else
            {
                _overlapAgents.push_back(  AgentPtr( agent ) );
            }
        }
#endif
    }

    void SpacePartition::sendGhostAgents( Rectangle<int> area, int dst )
    {
        std::cout << _id << " SEND to " << dst << "\n";
    }

    void SpacePartition::reciveGhostAgents( int src )
    {
        std::cout << _id << " RECV from " << src << "\n";
    }

    void SpacePartition::reduceGhostAgents( )
    {
        AgentsVector lateralAgents;
        AgentsVector verticalAgents;
        Overlap_st *lat = _mpiOverlap->getTopRight( );
        Overlap_st *ver = _mpiOverlap->getBottom( );
        for ( AgentsList::iterator it=_world->beginAgents( ); it!=_world->endAgents( ); it++ )+
        -

        // Send Right / Recv left
        if ( _mpiOverlap->isEven( ) )
        {
            if ( ( b = _mpiOverlap->getTopLeft( ) ) != 0 ) reciveGhostAgents( b->_n );
            if ( ( b = _mpiOverlap->getTopRight( ) ) != 0 ) sendGhostAgents( b->_local, b->_n );
            if ( ( b = _mpiOverlap->getTop( ) ) != 0 ) reciveGhostAgents( b->_n );
            if ( ( b = _mpiOverlap->getBottom( ) ) != 0 ) sendGhostAgents( b->_local, b->_n );
        }
        else
        {
            if ( ( b = _mpiOverlap->getTopRight( ) ) != 0 ) sendGhostAgents( b->_local, b->_n );
            if ( ( b = _mpiOverlap->getTopLeft( ) ) != 0 ) reciveGhostAgents( b->_n );
            if ( ( b = _mpiOverlap->getBottom( ) ) != 0 ) sendGhostAgents( b->_local, b->_n );
            if ( ( b = _mpiOverlap->getTop( ) ) != 0 ) reciveGhostAgents( b->_n );
        }
        abort( );
    }


#ifdef RGT
    void SpacePartition::reduceGhostAgents( )
    {
        int nBounds = _mpiOverlap->getNumOfBounds( );
        std::vector<Rectangle<int>> areas =  _mpiOverlap->getAreas( );
        AgentsVector ghostAgents[nBounds];
        for ( AgentsList::iterator it=_world->beginAgents( ); it!=_world->endAgents( ); it++ )
        {
            for ( size_t i= 0; i < nBounds; i++ )
            {
                if ( areas[i].contains( (*it)->getPosition( ) ) )
                {
                    ghostAgents[i].push_back( *it );
                    break;
                }
            }
        }
        for ( auto vec: ghostAgents )
            std::cout << _id << " -- " << vec.size() << "\n";
        

        abort();

        for ( int i= 0; i< 4; i++ )
        {
            if ( _nearby[i] != -1 )
            {
                for ( MpiFactory::TypesMap::iterator itType=MpiFactory::instance( )->beginTypes( ); itType!=MpiFactory::instance( )->endTypes( ); itType++ )
                {
                    MPI_Datatype * agentType = itType->second;
                    AgentsList list;
                    for ( auto agent: ghostAgents[i] )
                        if ( agent->isType( itType->first ) )
                            list.push_back( agent );
                    int sAgents = list.size( );
                    int rAgents;
                    MPI_Status status;
                    MPI_Sendrecv( &sAgents, 1, MPI_INTEGER, _nearby[i], eGhostAgent,
                                  &rAgents, 1, MPI_INTEGER, _nearby[i], eGhostAgent, MPI_COMM_WORLD, &status );
                    std::cout << _id << " -- " << _nearby[i] << " S/R " << itType->first << " " << sAgents << "/" << rAgents << "\n";
                    if (_id< _nearby[i])
                    {
                        sendAgentsList( list, _nearby[i], itType->second );
                        receiveAgentsList( rAgents, _nearby[i], itType->first, itType->second );
                    }
                    else
                    {
                        receiveAgentsList( rAgents, _nearby[i], itType->first, itType->second );
                        sendAgentsList( list, _nearby[i], itType->second );
                    }
                }
            }
        }
        abort( );
        for ( auto comm : _OverlapAreas.getNeighbors( ) )
        {
            int mpi_id = std::get<0>( comm );
            Rectangle<int> area = std::get<2>( comm ).intersection( _boundaries ); // Global coordinates
            for ( MpiFactory::TypesMap::iterator itType=MpiFactory::instance( )->beginTypes( ); itType!=MpiFactory::instance( )->endTypes( ); itType++ )
            {
                MPI_Datatype * agentType = itType->second;
                AgentsList list;
                for ( AgentsList::iterator it=_world->beginAgents( ); it!=_world->endAgents( ); it++ )
                {
                    AgentPtr agent = *it;
                    if ( area.contains( agent->getPosition( ) ) && agent->isType( itType->first ) )
                    {
                        if ( willBeRemoved( agent->getId( ) ) )
                        { std::cout << "What to do with agents that will be removed\n"; exit(0); }
                        list.push_back( agent );
                    }
                }
                for ( AgentsList::iterator it=_overlapAgents.begin( ); it!=_overlapAgents.end( ); it++ )
                {
                    AgentPtr agent = *it;
                    if ( area.contains( agent->getPosition( ) ) && agent->isType( itType->first ) ) 
                    {
                        if ( willBeRemoved( agent->getId( ) ) )
                        { std::cout << "What to do with agents that will be removed\n"; exit(0); }
                        list.push_back( agent );
                    }
                }
                int sAgents = list.size( );
                int rAgents;
                MPI_Status status;
                MPI_Sendrecv( &sAgents, 1, MPI_INTEGER, mpi_id, eGhostAgent,
                              &rAgents, 1, MPI_INTEGER, mpi_id, eGhostAgent, MPI_COMM_WORLD, &status );
                if (_id< mpi_id)
                {
                    sendAgentsList( list, mpi_id, itType->second );
                    receiveAgentsList( rAgents, mpi_id, itType->first, itType->second );
                }
                else
                {
                    receiveAgentsList( rAgents, mpi_id, itType->first, itType->second );
                    sendAgentsList( list, mpi_id, itType->second );
                }
            }
        }
    }
#endif

    void SpacePartition::receiveGhostAgents( const int & sectionIndex )
    {
        std::cout << "SpacePartition::receiveGhostAgents\n"; exit(1);
#ifdef RGT
        std::stringstream logName;
        logName << "MPI_agents_world_" << _id;
        log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " receive ghost agents for section index: " << sectionIndex );


        // we need to calculate how many neighbors will send data to this id
        std::vector<int> neighborsToUpdate;
        for ( size_t i=0; i<_neighbors.size( ); i++ )
        {
            if ( needsToReceiveData( _neighbors[i], sectionIndex ))
            {
                neighborsToUpdate.push_back( _neighbors[i] );
            }
        }

        // for each type of agent we will send the collection of agents of the particular type to neighbors
        for ( MpiFactory::TypesMap::iterator itType=MpiFactory::instance( )->beginTypes( ); itType!=MpiFactory::instance( )->endTypes( ); itType++ )
        {
            MPI_Datatype * agentType = itType->second;

            for ( size_t i=0; i<neighborsToUpdate.size( ); i++ )
            {
                log_DEBUG( logName.str( ), getWallTime( ) << " AA neigh to update: " << neighborsToUpdate[i] << " type: " << itType->first );
                AgentsList newGhostAgents;
                int numAgentsToReceive;
                MPI_Status status;

                int error = MPI_Recv( &numAgentsToReceive, 1, MPI_INTEGER, neighborsToUpdate[i], eNumGhostAgents, MPI_COMM_WORLD, &status );
                if ( error!=MPI_SUCCESS )
                {
                    std::stringstream oss;
                    oss << "SpacePartition::receiveGhostAgents - " << _id << " error in MPI_Recv: " << error;
                    throw Exception( oss.str( ) );
                }
                log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " has received message from " << neighborsToUpdate[i] << ", num ghost agents: " << numAgentsToReceive );
                for ( int j=0; j<numAgentsToReceive; j++ )
                {
                    void * package = MpiFactory::instance( )->createDefaultPackage( itType->first );
                    error = MPI_Recv( package, 1, *agentType, neighborsToUpdate[i], eGhostAgent, MPI_COMM_WORLD, &status );
                    if ( error!=MPI_SUCCESS )
                    {
                        std::stringstream oss;
                        oss << "SpacePartition::receiveGhostAgents - " << _id << " error in MPI_Recv: " << error;
                        throw Exception( oss.str( ) );
                    }
                    Agent * agent = MpiFactory::instance( )->createAndFillAgent( itType->first, package );
                    log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " has received ghost agent: " << agent << " number: " << j << " from: " << neighborsToUpdate[i] << " in section index: " << sectionIndex << " and step: " << _world->getCurrentStep( ) );
                    delete package;
                    agent->receiveVectorAttributes( neighborsToUpdate[i] );

                    // we must check if it is an update of an agent, or a ghost agent
                    bool worldOwnsAgent = false;
                    AgentsList::iterator it = getOwnedAgentPtr( agent->getId( ) );
                    if ( it!=_world->endAgents( ) )
                    {
                        log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " has received update of own agent: " << *it << " in step: " << _world->getCurrentStep( ) );
                        _world->eraseAgent( it );
                        _world->addAgent( agent, false );
                        worldOwnsAgent = true;
                    }
                    if ( !worldOwnsAgent )
                    {
                        newGhostAgents.push_back( std::shared_ptr<Agent>( agent ));
                    }
                }
                log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " num ghost agents sent for neighbor: " << neighborsToUpdate[i] << " in section index: " << sectionIndex << ": " << newGhostAgents.size( ) );
                // if the agent is in the zone to be updated, remove it
                Rectangle<int> overlapZone = getOverlap( neighborsToUpdate[i], sectionIndex );
                AgentsList::iterator it=_overlapAgents.begin( );
                while( it!=_overlapAgents.end( ) )
                {
                    Agent * agent = it->get( );
                    if ( agent->isType( itType->first ))
                    {
                        // si l'agent no està en zona que s'ha d'actualitzar, continuar
                        if ( overlapZone.contains( (*it )->getPosition( )-_boundaries._origin ))
                        {
                            log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " in section index: " << sectionIndex << " with overlap zone: " << overlapZone << " erasing agent: " << *it );
                            it = _overlapAgents.erase( it );
                        }
                        else
                        {
                            log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " in section index: " << sectionIndex <<  " with overlap zone: " << overlapZone << " maintaining agent: " << *it );
                            it++;
                        }
                    }
                    else
                    {
                        it++;
                    }
                }
                // afterwards we will add the new ghost agents
                for ( it=newGhostAgents.begin( ); it!=newGhostAgents.end( ); it++ )
                {
                    AgentPtr agent = *it;
                    log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " in section index: " << sectionIndex << " adding ghost agent: " << agent );
                    _overlapAgents.push_back( agent );
                }
            }
        }
        log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " receive ghost agents for section index: " << sectionIndex << " finished" );
#endif
    }

    void SpacePartition::receiveAgents( const int & sectionIndex )
    {
        std::cout << "SpacePartition::receiveAgents\n"; exit(1);
#ifdef RGT
        std::stringstream logName;
        logName << "MPI_agents_world_" << _id;
        log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " receiving agents for section index: " << sectionIndex );
        for ( MpiFactory::TypesMap::iterator itType=MpiFactory::instance( )->beginTypes( ); itType!=MpiFactory::instance( )->endTypes( ); itType++ )
        {
            log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " receiveAgent - checking mpi type: " << itType->first );
            MPI_Datatype * agentType = itType->second;

            for ( size_t i=0; i<_neighbors.size( ); i++ )
            {
                log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " receiveAgent - checking mpi type: " << itType->first << " for neighbor: " << _neighbors[i] );

                int numAgentsToReceive;
                MPI_Status status;
                int error = MPI_Recv( &numAgentsToReceive, 1, MPI_INTEGER, _neighbors[i], eNumAgents, MPI_COMM_WORLD, &status );
                if ( error!=MPI_SUCCESS )
                {
                    std::stringstream oss;
                    oss << "SpacePartition::receiveAgents - " << _id << " error in MPI_Recv: " << error;
                    throw Exception( oss.str( ) );
                }
                log_DEBUG( logName.str( ), getWallTime( ) <<  " receiveAgents - received message from " << _neighbors[i] << ", num agents: " << numAgentsToReceive );
                for ( int j=0; j<numAgentsToReceive; j++ )
                {
                    void * package = MpiFactory::instance( )->createDefaultPackage( itType->first );
                    error = MPI_Recv( package, 1, *agentType, _neighbors[i], eAgent, MPI_COMM_WORLD, &status );
                    if ( error!=MPI_SUCCESS )
                    {
                        std::stringstream oss;
                        oss << "SpacePartition::receiveAgents - " << _id << " error in MPI_Recv: " << error;
                        throw Exception( oss.str( ) );
                    }
                    Agent * agent = MpiFactory::instance( )->createAndFillAgent( itType->first, package );
                    delete package;
                    agent->receiveVectorAttributes( _neighbors[i] );
                    _world->addAgent( agent, true );
                    log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " receiveAgents - received agent: " << agent << " number: " << j << " from: " << _neighbors[i] );
                }
            }
        }
#endif
    }

    void SpacePartition::receiveOverlapData( const int & sectionIndex, const bool & entireOverlap )
    {
        std::cout << "SpacePartition::receiveOverlapData\n"; exit(1);
#ifdef RGT
        std::stringstream logName;
        logName << "MPI_raster_world_" << _id;
        log_DEBUG( logName.str( ), getWallTime( ) << " step: "  << "/" << sectionIndex << _world->getCurrentStep( ) << " receiveOverlapData" );

        // we need to calculate how many neighbors will send data to this id
        std::vector<int> neighborsToUpdate;
        for ( size_t i=0; i<_neighbors.size( ); i++ )
        {
            if ( needsToReceiveData( _neighbors[i], sectionIndex ))
            {
                neighborsToUpdate.push_back( _neighbors[i] );
            }
        }

        // for each raster, we receive data from all the active neighbors
        for ( size_t d=0; d<_world->getNumberOfRasters( ); d++ )
        {
            if ( !_world->rasterExists( d ) || !_world->isRasterDynamic( d ))
            {
                continue;
            }

            log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << "/" << sectionIndex << " receiving raster: " << d );
            for ( size_t i=0; i<neighborsToUpdate.size( ); i++ )
            {
                MpiOverlap* receive = new MpiOverlap;
                // TODO move to index
                receive->_rasterName = _world->getRasterName( d );
                if ( entireOverlap )
                {
                    receive->_overlap = getOverlap( neighborsToUpdate[i], sectionIndex );
                    log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << "/" << sectionIndex << " receive entire overlap: " << receive->_overlap << " from " << neighborsToUpdate[i] );
                }
                else
                {
                    receive->_overlap = getExternalOverlap( neighborsToUpdate[i] );
                    log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << "/" << sectionIndex << " receive external overlap: " << receive->_overlap << " from " << neighborsToUpdate[i] );
                }
                receive->_data.resize( receive->_overlap._size._width*receive->_overlap._size._height );
                log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << "/" << sectionIndex << " will receive overlap from: " << neighborsToUpdate[i] << " with size: " << receive->_data.size( ) << " and zone: " << receive->_overlap );
                MPI_Irecv( &receive->_data[0], receive->_data.size( ), MPI_INTEGER, neighborsToUpdate[i], eRasterData, MPI_COMM_WORLD, &receive->_request );
                _receiveRequests.push_back( receive );
            }
            log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << "/" << sectionIndex << " raster: " << d << " received" );
        }
        log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << "/" << sectionIndex << " receiveOverlapData ended" );
#endif
    }

    void SpacePartition::receiveMaxOverlapData( )
    {
        std::cout << "SpacePartition::receiveMaxOverlapData\n"; exit(1);
#ifdef RGT
        std::stringstream logName;
        logName << "MPI_raster_world_" << _id;
        log_DEBUG( logName.str( ), getWallTime( ) << " step: "  << _world->getCurrentStep( ) << " receiveMaxOverlapData" );
        // for each raster, we receive data from all the active neighbors
        for ( size_t d=0; d<_world->getNumberOfRasters( ); d++ )
        {
            if ( !_world->rasterExists( d ) || !_world->isRasterDynamic( d ))
            {
                continue;
            }
            log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " receiving max raster: " << d );
            for ( size_t i=0; i<_neighbors.size( ); i++ )
            {
                MpiOverlap* receive = new MpiOverlap;
                receive->_rasterName = _world->getRasterName( d );
                receive->_overlap = getExternalOverlap( _neighbors[i] );
                receive->_data.resize( receive->_overlap._size._width*receive->_overlap._size._height );

                log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " will receive max overlap to: " << _neighbors[i] << " with size: " << receive->_data.size( ) << " and zone: " << receive->_overlap << " from " << _neighbors[i] );
                MPI_Irecv( &receive->_data[0], receive->_data.size( ), MPI_INTEGER, _neighbors[i], eRasterMaxData, MPI_COMM_WORLD, &receive->_request );
                _receiveRequests.push_back( receive );
            }
            log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( )  << " raster: " << d << " max data received" );
        }
        log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " receiveMaxOverlapData ended" );
#endif
    }

    void SpacePartition::clearRequests( bool updateMaxValues )
    {
        std::cout << "SpacePartition::clearRequests\n"; exit(1);
        std::stringstream logName;
        logName << "MPI_raster_world_" << _id;
        int finished = 0;

        // MPI_Isend
        while( _sendRequests.size( )!=0 )
        {
            std::list<MpiOverlap*>::iterator it=_sendRequests.begin( );
            while( it!=_sendRequests.end( ) )
            {
                MpiOverlap* send = *it;
                MPI_Status status;
                MPI_Test( &send->_request, &finished, &status );
                if ( finished )
                {
                    it = _sendRequests.erase( it );
                    log_DEBUG( logName.str( ), getWallTime( ) << " request finished, lacking: " << _sendRequests.size( ) << " requests" );
                    delete send;
                }
            }
        }

        // MPI_Irecv
        while( _receiveRequests.size( )!=0 )
        {
            std::list<MpiOverlap*>::iterator it=_receiveRequests.begin( );
            while( it!=_receiveRequests.end( ) )
            {
                MpiOverlap* receive = *it;
                MPI_Status status;
                MPI_Test( &receive->_request, &finished, &status );
                if ( finished )
                {
                    it = _receiveRequests.erase( it );
                    log_DEBUG( logName.str( ), getWallTime( ) << " receive request finished, lacking: " << _receiveRequests.size( ) << " requests" );

                    const Rectangle<int> & overlapZone = receive->_overlap;
                    log_DEBUG( logName.str( ), getWallTime( ) << " processing request for overlap: " << overlapZone << " raster: " << receive->_rasterName << " for max values" );
                    for ( size_t i=0; i<receive->_data.size( ); i++ )
                    {
                        Point2D<int> index( overlapZone._origin._x+i%overlapZone._size._width, overlapZone._origin._y+i/overlapZone._size._width );
                        if ( updateMaxValues )
                        {
                            log_EDEBUG( logName.str( ), "\t" << getWallTime( ) << " step: " << _world->getCurrentStep( ) << " receive index: " << index << " max value: " << receive->_data.at( i ));
                            _world->getDynamicRaster( receive->_rasterName ).setMaxValue( index, receive->_data.at( i ));
                        }
                        else
                        {
                            log_EDEBUG( logName.str( ), "\t" << getWallTime( ) << " step: " << _world->getCurrentStep( ) << " receive index: " << index << " current value: " << receive->_data.at( i ));
                            _world->getDynamicRaster( receive->_rasterName ).setValue( index, receive->_data.at( i ));
                        }
                    }

                    if ( updateMaxValues )
                    {
                        log_DEBUG( logName.str( ), getWallTime( ) << " receive request finished for max values of raster: " <<  receive->_rasterName );
                    }
                    else
                    {
                        log_DEBUG( logName.str( ), getWallTime( ) << " receive request finished for current values of raster: " <<  receive->_rasterName );
                    }
                    delete receive;
                }
            }
        }

    }

    int SpacePartition::getIdFromPosition( const Point2D<int> & position )
    {
        std::cout << "SpacePartition::getIdFromPosition\n"; exit(1);
#ifdef KKK
        Point2D<int> nodePosition( position._x/_ownedArea._size._width, position._y/_ownedArea._size._height );
        return nodePosition._y*sqrt( _numTasks )+nodePosition._x;
#endif
    }

    Point2D<int> SpacePartition::getPositionFromId( const int & id ) const
    {
        int worldsPerRow = sqrt( _numTasks );
        Point2D<int> worldPos( id%worldsPerRow, id/worldsPerRow );
        return worldPos;
    }

    int SpacePartition::getNeighborIndex( const int & id )
    {
        std::cout << "SpacePartition::getNeighborIndex\n"; exit(1);
#ifdef RGT
        for ( size_t i=0; i<_neighbors.size( ); i++ )
        {
            if ( _neighbors[i] == id )
            {
                return i;
            }
        }
        std::stringstream oss;
        oss << "SpacePartition::getNeighborIndex- " << _id << " looking for neighbor with id: " << id << " not found";
        throw Exception( oss.str( ) );
#endif
    }

    void SpacePartition::executeAgents( )
    {
        std::cout << _id << " n_a=" << _world->getNumberOfAgents( ) << " ov_a=" << _overlapAgents.size( ) << "\n";
#ifndef ORIG
#ifdef RGT
        // Reset Executed Agents map
        _executedAgentsHash.clear( );

        AgentsVector agentsToExecute[4];
        for ( AgentsList::iterator it=_world->beginAgents( ); it!=_world->endAgents( ); it++ )
            if ( ! hasBeenExecuted( (*it)->getId( ) ) )
                agentsToExecute[_OverlapAreas.getSection( (*it)->getPosition( ) )].push_back( *it );


        for ( int sectionIndex=0; sectionIndex<4; sectionIndex++ )
        {
            stepSection( sectionIndex, agentsToExecute[sectionIndex] );
        
            agentsToExecute[sectionIndex].clear( );
        }
        std::cout << _id << " SpacePartition::executeAgents\n"; MPI_Barrier( MPI_COMM_WORLD ); MPI_Finalize( ); exit(0);
#endif

#else
        for ( int sectionIndex=0; sectionIndex<4; sectionIndex++ )
        {
            // section index doesn't matter if is the entire overlap
            // TODO refactor? we are sending 4 times all the info
            sendOverlapZones( sectionIndex, false );
            receiveOverlapData( sectionIndex, false );
        }

        std::stringstream logName;
        logName << "simulation_" << getId( );

        log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " has executed update overlap" );
        _executedAgentsHash.clear( );

        std::stringstream logNameMpi;
        logNameMpi << "simulation_" << _id;

        log_DEBUG( logName.str( ), getWallTime( ) << " step: " << _world->getCurrentStep( ) << " executing sections" );
        for ( int sectionIndex=0; sectionIndex<4; sectionIndex++ )
        {
            stepSection( sectionIndex );
            log_DEBUG( logNameMpi.str( ), getWallTime( ) << " executing step: " << _world->getCurrentStep( ) << " and section: " << sectionIndex << " has been executed" );

            receiveAgents( sectionIndex );
            log_DEBUG( logNameMpi.str( ), getWallTime( ) << " executing step: " << _world->getCurrentStep( ) << " and section: " << sectionIndex << " has received agents" );

            sendGhostAgents( sectionIndex );
            log_DEBUG( logNameMpi.str( ), getWallTime( ) << " executing step: " << _world->getCurrentStep( ) << " and section: " << sectionIndex << " sent ghosts" );

            receiveGhostAgents( sectionIndex );
            log_DEBUG( logNameMpi.str( ), getWallTime( ) << " executing step: " << _world->getCurrentStep( ) << " and section: " << sectionIndex << " received ghosts" );

            MPI_Barrier( MPI_COMM_WORLD );

            sendOverlapZones( sectionIndex );
            log_DEBUG( logNameMpi.str( ), getWallTime( ) << " executing step: " << _world->getCurrentStep( ) << " and section: " << sectionIndex << " sent overlap" );

            receiveOverlapData( sectionIndex );
            log_DEBUG( logNameMpi.str( ), getWallTime( ) << " executing step: " << _world->getCurrentStep( ) << " and section: " << sectionIndex << " received overlap" );

            MPI_Barrier( MPI_COMM_WORLD );

            clearRequests( false );
        }
#endif
    }

    void SpacePartition::initOverlappingData( )
    {
        // we need to send the agents and data in overlap zone to adjacent computer nodes
#ifdef RGT
        sendMaxOverlapZones( );
        std::cout << "SpacePartition::initOverlappingData 2\n"; exit(1);
        receiveMaxOverlapData( );
        // all nodes must finish receiving max values before receiving current values
        MPI_Barrier( MPI_COMM_WORLD );
        clearRequests( true );
        for ( int sectionIndex=0; sectionIndex<4; sectionIndex++ )
        {

            // section index doesn't matter if is the entire overlap
            // TODO refactor? we are sending 4 times all the info
            sendOverlapZones( sectionIndex, false );
            receiveOverlapData( sectionIndex, false );
            // all nodes must finish receiving max values before receiving current values
            MPI_Barrier( MPI_COMM_WORLD );
            clearRequests( false );


            sendGhostAgents( sectionIndex );
            receiveGhostAgents( sectionIndex );
        }
#else
        for ( size_t d=0; d<_world->getNumberOfRasters( ); d++ )
            if ( _world->rasterExists( d )  && _world->isRasterDynamic( d ) )
                _rasters.push_back( &_world->getDynamicRaster( d ) );
        reduceOverlapZones( );
        reduceGhostAgents( );
#endif


    }

    void SpacePartition::abort( )
    {
        MPI_Barrier( MPI_COMM_WORLD );
        MPI_Finalize( );
        exit(0);
    }

    void SpacePartition::finish( )
    {
        std::stringstream logName;
        logName << "simulation_" << _id;

        _serializer.finish( );

        log_INFO( logName.str( ), getWallTime( ) << " simulation finished" );
        if ( _finalize )
        {
            MpiFactory::instance( )->cleanTypes( );
            MPI_Finalize( );
        }
    }

    const Rectangle<int> & SpacePartition::getBoundaries( ) const
    {
        return _boundaries;
    }

    bool SpacePartition::isCorner( const int & neighbor ) const
    {
        std::cout << "SpacePartition::isCorner\n"; exit(1);
        Point2D<int> worldPos = getPositionFromId( neighbor );
        Point2D<int> diff = worldPos - _worldPos;
        if ( std::abs( diff._x )==1 && std::abs( diff._y )==1 )
        {
            return true;
        }
        return false;
    }

    Rectangle<int> SpacePartition::getInternalOverlap( const int & id ) const
    {
        std::cout << "SpacePartition::getInternalOverlap\n"; exit(1);
        Point2D<int> diff = getPositionFromId( id )-_worldPos;
        // left

        Rectangle<int> result;

        // origin
        if ( diff._x==-1 )
        {
            if ( _ownedArea._origin._x==0 )
            {
                result._origin._x = 0;
            }
            else
            {
                result._origin._x = _overlap;
            }
        }
        else if ( diff._x==0 )
        {
            if ( _ownedArea._origin._x==0 )
            {
                result._origin._x = 0;
            }
            else
            {
                result._origin._x = _overlap;
            }
        }
        else
        {
            // if left border just remove an overlap
            if ( _ownedArea._origin._x==0 )
            {
                result._origin._x = _ownedArea._size._width - _overlap;
            }
            // else sum and remove an overlap
            else
            {
                result._origin._x = _ownedArea._size._width;
            }
        }

        if ( diff._y==-1 )
        {
            if ( _ownedArea._origin._y==0 )
            {
                result._origin._y = 0;
            }
            else
            {
                result._origin._y = _overlap;
            }
        }
        else if ( diff._y==0 )
        {
            if ( _ownedArea._origin._y==0 )
            {
                result._origin._y = 0;
            }
            else
            {
                result._origin._y = _overlap;
            }
        }
        else
        {
            // if top border just remove an overlap
            if ( _ownedArea._origin._y==0 )
            {
                result._origin._y = _ownedArea._size._height - _overlap;
            }
            // else sum and remove an overlap
            else
            {
                result._origin._y = _ownedArea._size._height;
            }
        }

        // size
        if ( diff._x!=0 )
        {
            result._size._width = _overlap;
        }
        else
        {
            result._size._width = _ownedArea._size._width;
        }

        if ( diff._y!=0 )
        {
            result._size._height = _overlap;
        }
        else
        {
            result._size._height = _ownedArea._size._height;
        }
        std::stringstream logName;
        logName << "MPI_raster_world_" << _id;
        log_DEBUG( logName.str( ), getWallTime( ) << " internal overlap with: " << id << " and diff: " << diff << " - " << result );
        return result;
    }

    Rectangle<int> SpacePartition::getExternalOverlap( const int & id ) const
    {
        std::cout << "SpacePartition::getExternalOverlap\n"; exit(1);
        Point2D<int> diff = getPositionFromId( id )-_worldPos;
        // left

        Rectangle<int> result;

        // origin
        if ( diff._x==-1 )
        {
            result._origin._x = 0;
        }
        else if ( diff._x==0 )
        {
            if ( _ownedArea._origin._x == 0 )
            {
                result._origin._x = 0;
            }
            else
            {
                result._origin._x = _overlap;
            }
        }
        else
        {
            // if left border it doesn't have a left overlap
            if ( _ownedArea._origin._x==0 )
            {
                result._origin._x = _ownedArea._size._width;
            }
            // else sum an overlap
            else
            {
                result._origin._x = _ownedArea._size._width + _overlap;
            }
        }

        if ( diff._y==-1 )
        {
            result._origin._y = 0;
        }
        else if ( diff._y==0 )
        {
            if ( _ownedArea._origin._y == 0 )
            {
                result._origin._y = 0;
            }
            else
            {
                result._origin._y = _overlap;
            }
        }
        else
        {
            // if lower border it doesn't have a lower overlap
            if ( _ownedArea._origin._y==0 )
            {
                result._origin._y = _ownedArea._size._height;
            }
            // else sum an overlap
            else
            {
                result._origin._y = _ownedArea._size._height + _overlap;
            }
        }

        // size
        if ( diff._x!=0 )
        {
            result._size._width = _overlap;
        }
        else
        {
            result._size._width = _ownedArea._size._width;
        }

        if ( diff._y!=0 )
        {
            result._size._height = _overlap;
        }
        else
        {
            result._size._height = _ownedArea._size._height;
        }

        std::stringstream logName;
        logName << "MPI_raster_world_" << _id;
        log_DEBUG( logName.str( ), getWallTime( ) << " external overlap with: " << id << " and diff: " << diff << " - " << result );
        return result;
    }

    Rectangle<int> SpacePartition::getOverlap( const int & id, const int & sectionIndex ) const
    {
        std::cout << "SpacePartition::getOverlap 2\n"; exit(1);
        Point2D<int> diff = getPositionFromId( id )-_worldPos;
        // left

        Rectangle<int> result;

        // origin
        if ( diff._x==-1 )
        {
            result._origin._x = 0;
        }
        else if ( diff._x==0 )
        {
            if ( sectionIndex==0 || sectionIndex==2 )
            {
                result._origin._x = 0;
            }
            else
            {
                result._origin._x = _ownedArea._size._width/2;
                if ( _ownedArea._origin._x+_ownedArea._size._width!=_world->getConfig( ).getSize( )._width )
                {
                    result._origin._x -= _overlap;
                }
            }
        }
        else
        {
            // if left border just remove an overlap
            if ( _ownedArea._origin._x==0 )
            {
                result._origin._x = _ownedArea._size._width - _overlap;
            }
            // else sum and remove an overlap
            else
            {
                result._origin._x = _ownedArea._size._width;
            }
        }

        if ( diff._y==-1 )
        {
            result._origin._y = 0;
        }
        else if ( diff._y==0 )
        {
            if ( sectionIndex==0 || sectionIndex==1 )
            {
                result._origin._y = 0;
            }
            else
            {
                result._origin._y = _ownedArea._size._height/2;
                if ( _ownedArea._origin._y+_ownedArea._size._height!=_world->getConfig( ).getSize( )._height )
                {
                    result._origin._y -= _overlap;
                }
            }
        }
        else
        {
            // if top border just remove an overlap
            if ( _ownedArea._origin._y==0 )
            {
                result._origin._y = _ownedArea._size._height - _overlap;
            }
            // else sum and remove an overlap
            else
            {
                result._origin._y = _ownedArea._size._height;
            }
        }

        // size
        if ( diff._x!=0 )
        {
            result._size._width = _overlap*2;
        }
        else
        {
            // this result is = to boundaries + 1 overlap
            result._size._width = _ownedArea._size._width/2;
            // border
            if ( _ownedArea._origin._x!=0 )
            {
                result._size._width += _overlap;
            }
            else
            {
                if ( sectionIndex!=0 && sectionIndex!=2 )
                {
                    result._size._width += _overlap;
                }
            }

            if ( _ownedArea._origin._x+_ownedArea._size._width!=_world->getConfig( ).getSize( )._width )
            {
                result._size._width += _overlap;
            }
            else
            {
                if ( sectionIndex!=1 && sectionIndex!=3 )
                {
                    result._size._width += _overlap;
                }
            }
        }

        if ( diff._y!=0 )
        {
            result._size._height = _overlap*2;
        }
        else
        {
            result._size._height = _ownedArea._size._height/2;
            if ( _ownedArea._origin._y!=0 )
            {
                result._size._height += _overlap;
            }
            else
            {
                if ( sectionIndex!=0 && sectionIndex!=1 )
                {
                    result._size._height += _overlap;
                }
            }

            if ( _ownedArea._origin._y+_ownedArea._size._height!=_world->getConfig( ).getSize( )._height )
            {
                result._size._height += _overlap;
            }
            else
            {
                if ( sectionIndex!=2 && sectionIndex!=3 )
                {
                    result._size._height += _overlap;
                }
            }
        }
        std::stringstream logName;
        logName << "MPI_raster_world_" << _id;
        log_DEBUG( logName.str( ), getWallTime( ) << " overlap with: " << id << " and diff: " << diff << " in section Index: " << sectionIndex << " - " << result );
        return result;
    }

    bool SpacePartition::needsToBeUpdated( const int & id, const int & sectionIndex )
    {
        std::cout << "SpacePartition::needsToBeUpdated\n"; exit(1);
        Point2D<int> diff = getPositionFromId( id )-_worldPos;
        switch( sectionIndex )
        {
            case 0:
                // top left
                if ( diff._x==-1 && diff._y==-1 )
                {
                    return true;
                }
                // left
                if ( diff._x==-1 && diff._y==0 )
                {
                    return true;
                }
                // top
                if ( diff._x==0 && diff._y==-1 )
                {
                    return true;
                }
                return false;
                break;
            case 1:
                // top
                if ( diff._x==0 && diff._y==-1 )
                {
                    return true;
                }
                // top right
                if ( diff._x==1 && diff._y==-1 )
                {
                    return true;
                }
                // right
                if ( diff._x==1 && diff._y==0 )
                {
                    return true;
                }
                return false;
                break;
            case 2:
                // bottom
                if ( diff._x==0 && diff._y==1 )
                {
                    return true;
                }
                // bottom left
                if ( diff._x==-1 && diff._y==1 )
                {
                    return true;
                }
                // left
                if ( diff._x==-1 && diff._y==0 )
                {
                    return true;
                }
                return false;
                break;
            case 3:
                // right
                if ( diff._x==1 && diff._y==0 )
                {
                    return true;
                }
                // bottom right
                if ( diff._x==1 && diff._y==1 )
                {
                    return true;
                }
                // bottom
                if ( diff._x==0 && diff._y==1 )
                {
                    return true;
                }
                return false;
                break;
            default:
                std::stringstream oss;
                oss << "SpacePartition::needsToBeUpdated - wrong section index: " << sectionIndex;
                throw Engine::Exception( oss.str( ) );
        }
        return false;
    }

    bool SpacePartition::needsToReceiveData( const int & id, const int & sectionIndex )
    {
        std::cout << "SpacePartition::needsToReceiveData\n"; exit(1);
        Point2D<int> diff = getPositionFromId( id )-_worldPos;
        switch( sectionIndex )
        {
            case 0:
                // bottom
                if ( diff._x==0 && diff._y==1 )
                {
                    return true;
                }
                // right
                if ( diff._x==1 && diff._y==0 )
                {
                    return true;
                }
                // bottom right
                if ( diff._x==1 && diff._y==1 )
                {
                    return true;
                }
                return false;
                break;
            case 1:
                // left
                if ( diff._x==-1 && diff._y==0 )
                {
                    return true;
                }
                // bottom left
                if ( diff._x==-1 && diff._y==1 )
                {
                    return true;
                }
                // bottom
                if ( diff._x==0 && diff._y==1 )
                {
                    return true;
                }
                return false;
                break;
            case 2:
                // top
                if ( diff._x==0 && diff._y==-1 )
                {
                    return true;
                }
                // top right
                if ( diff._x==1 && diff._y==-1 )
                {
                    return true;
                }
                // right
                if ( diff._x==1 && diff._y==0 )
                {
                    return true;
                }
                return false;
                break;
            case 3:
                // top left
                if ( diff._x==-1 && diff._y==-1 )
                {
                    return true;
                }
                // left
                if ( diff._x==-1 && diff._y==0 )
                {
                    return true;
                }
                // top
                if ( diff._x==0 && diff._y==-1 )
                {
                    return true;
                }
                return false;
                break;
            default:
                std::stringstream oss;
                oss << "SpacePartition::needsToReceiveData - wrong section index: " << sectionIndex;
                throw Engine::Exception( oss.str( ) );
        }
        return false;
    }

    const int & SpacePartition::getOverlap( ) const
    {
        return _overlap;
    }

    bool SpacePartition::hasBeenExecuted( const std::string & id ) const
    {
        if ( _executedAgentsHash.find( id ) == _executedAgentsHash.end( ) )
        {
            return false;
        }
        return true;
    }

    void SpacePartition::agentAdded( AgentPtr agent, bool executedAgent )
    {
        _executedAgentsHash.insert( make_pair( agent->getId( ), agent ));
        AgentsList::iterator it = getGhostAgentPtr( agent->getId( ) );
        if ( it!=_overlapAgents.end( ) )
        {
            _overlapAgents.erase( it );
        }
        if ( !executedAgent )
        {
            return;
        }
    }

    AgentsList::iterator SpacePartition::getGhostAgentPtr( const std::string & id )
    {
        for ( AgentsList::iterator it=_overlapAgents.begin( ); it!=_overlapAgents.end( ); it++ )
        {
            if ( (*it )->getId( ).compare( id )==0 )
            {
                return it;
            }
        }
        return _overlapAgents.end( );
    }

    Agent *SpacePartition::getGhostAgent( const std::string & id )
    {
        for ( AgentsList::iterator it=_overlapAgents.begin( ); it!=_overlapAgents.end( ); it++ )
        {
            if ( (*it )->getId( ).compare( id )==0 )
            {
                return it->get( );
            }
        }
        return 0;
    }

    Agent * SpacePartition::getAgent( const std::string & id )
    {
        std::cout << "SpacePartition::getAgent\n"; exit(1);
        Agent *agent = getOwnedAgent( id );
        if ( agent != 0 && agent->exists( ) )
            return agent;
        agent = getGhostAgent( id );
        if ( agent != 0 && agent->exists( ) )
            return agent;
        return 0;
    }

    void SpacePartition::removeAgent( Agent * agent )
    {
        std::cout << "SpacePartition::removeAgent\n"; exit(1);
        AgentsList::iterator it = getOwnedAgentPtr( agent->getId( ) );
        if ( it==_world->endAgents( ) )
        {
            std::stringstream oss;
            oss << "SpacePartition::removeAgent - agent: " << agent << " not found";
            throw Exception( oss.str( ) );
        }
        // TODO it is not needed if it has modified position, as it is already done after the executed of a given agent step
        //sendDeleteOverlapAgent( it, agent->getPosition( ) );
        _removedAgents.push_back( *it );
    }

    void SpacePartition::removeAgents( )
    {
        std::cout << "SpacePartition::removeAgent\n"; exit(1);
        AgentsList::iterator it=_removedAgents.begin( );
        while( it!=_removedAgents.end( ) )
        {
            Agent * agent = it->get( );
            AgentsList::iterator itAg = getOwnedAgentPtr( agent->getId( ) );
            if ( itAg==_world->endAgents( ) )
            {
                std::stringstream oss;
                oss << "SpacePartition::removeAgents - agent: " << agent << " not found";
                throw Exception( oss.str( ) );
                return;
            }
            _world->eraseAgent( itAg );
            it = _removedAgents.erase( it );
        }
        _removedAgents.clear( );
    }

    AgentsVector SpacePartition::getAgent( const Point2D<int> & position, const std::string & type )
    {
        std::cout << "SpacePartition::getAgent\n"; exit(1);
        AgentsVector result;
        for ( AgentsList::iterator it=_world->beginAgents( ); it!=_world->endAgents( ); it++ )
        {
            AgentPtr agent = *it;
            if ( agent->getPosition( ).isEqual( position ))
            {
                if ( type.compare( "all" )==0 || agent->isType( type ))
                {
                    result.push_back( agent );
                }
            }
        }
        for ( AgentsList::iterator it=_overlapAgents.begin( ); it!=_overlapAgents.end( ); it++ )
        {
            AgentPtr agent = *it;
            if ( agent->getPosition( ).isEqual( position ))
            {
                if ( type.compare( "all" )==0 || agent->isType( type ))
                {
                    result.push_back( agent );
                }
            }
        }
        return result;
    }

    Agent *SpacePartition::getOwnedAgent( const std::string & id )
    {
        for ( AgentsList::iterator it=_world->beginAgents( ); it!=_world->endAgents( ); it++ )
        {
            if ( (*it )->getId( ).compare( id )==0 )
            {
                return it->get( );
            }
        }
        return 0;
    }

    AgentsList::iterator SpacePartition::getOwnedAgentPtr( const std::string & id )
    {
        for ( AgentsList::iterator it=_world->beginAgents( ); it!=_world->endAgents( ); it++ )
        {
            if ( (*it )->getId( ).compare( id )==0 )
            {
                return it;
            }
        }
        return _world->endAgents( );
    }

    bool SpacePartition::willBeRemoved( const std::string & id )
    {
        for ( AgentsList::iterator it=_removedAgents.begin( ); it!=_removedAgents.end( ); it++ )
        {
            if ( (*it )->getId( ).compare( id )==0 )
            {
                return true;
            }
        }
        return false;
    }

    const Rectangle<int> & SpacePartition::getOwnedArea( ) const
    {
        return _ownedArea;
    }

    Point2D<int> SpacePartition::getRealPosition( const Point2D<int> & globalPosition ) const
    {
        return globalPosition-_boundaries._origin;
    }

    Point2D<int> SpacePartition::getRandomPosition( ) const
    {
        Engine::Point2D<int> pos( GeneralState::statistics( ).getUniformDistValue( 0, _ownedArea._size._width-1 ),
                                  GeneralState::statistics( ).getUniformDistValue( 0, _ownedArea._size._height-1 ) );
        pos += _ownedArea._origin;
        return pos;
    }

    double SpacePartition::getWallTime( ) const
    {
        return MPI_Wtime( ) - _initialTime;
    }

    size_t SpacePartition::getNumberOfTypedAgents( const std::string & type ) const
    {
        size_t n = 0;
        AgentsList::iterator it = _world->beginAgents( );
        while( it!=_world->endAgents( ) )
        {
            AgentPtr agent = *it;
            if ( agent->isType( type ) ) n++;
            it++;
        }
        size_t N;
        MPI_Allreduce( &n, &N, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD );
        return N;
    }

    void SpacePartition::addStringAttribute( const std::string & type, const std::string & key, const std::string & value )
    {
        _serializer.addStringAttribute( type, key, value );
    }

    void SpacePartition::addIntAttribute( const std::string & type, const std::string & key, int value )
    {
        _serializer.addIntAttribute( type, key, value );
    }

    void SpacePartition::addFloatAttribute( const std::string & type, const std::string & key, float value )
    {
        _serializer.addFloatAttribute( type, key, value );
    }

    void SpacePartition::serializeAgents( const int & step )
    {
        _serializer.serializeAgents( step, _world->beginAgents( ), _world->endAgents( ) );
    }

    void SpacePartition::serializeRasters( const int & step )
    {
        _serializer.serializeRasters( step );
    }

    int SpacePartition::countNeighbours( Agent * target, const double & radius, const std::string & type )
    {
        std::cout << "SpacePartition::countNeighbours 1\n"; exit(1);
        int numAgents = for_each( _world->beginAgents( ), _world->endAgents( ), aggregatorCount<std::shared_ptr<Agent> >( radius, *target, type ))._count;
        int numOverlapAgents = for_each( _overlapAgents.begin( ), _overlapAgents.end( ), aggregatorCount<std::shared_ptr<Agent> >( radius, *target, type ))._count;
        return numAgents+numOverlapAgents;
    }

    AgentsVector SpacePartition::getNeighbours( Agent * target, const double & radius, const std::string & type )
    {
        AgentsVector agentsVector = for_each( _world->beginAgents( ), _world->endAgents( ), aggregatorGet<std::shared_ptr<Agent> >( radius, *target, type ))._neighbors;
        AgentsVector overlapAgentsVector =  for_each( _overlapAgents.begin( ), _overlapAgents.end( ), aggregatorGet<std::shared_ptr<Agent> >( radius, *target, type ))._neighbors;
        std::copy( overlapAgentsVector.begin( ), overlapAgentsVector.end( ), std::back_inserter( agentsVector ));
        std::random_shuffle( agentsVector.begin( ), agentsVector.end( ) );
        return agentsVector;
    }

    void SpacePartition::setValue( DynamicRaster & raster, const Point2D<int> & position, int value )
    {
        raster.setValue( getRealPosition( position ), value );
    }

    int SpacePartition::getValue( const DynamicRaster & raster, const Point2D<int> & position ) const
    {
        return raster.getValue( getRealPosition( position ) );
    }

    void SpacePartition::setMaxValue( DynamicRaster & raster, const Point2D<int> & position, int value )
    {
        std::cout << "SpacePartition::setMaxValue\n"; exit(1);
        raster.setMaxValue( getRealPosition( position ), value );
    }

    int SpacePartition::getMaxValue( const DynamicRaster & raster, const Point2D<int> & position ) const
    {
        std::cout << "SpacePartition::getMaxValue\n"; exit(1);
        return raster.getMaxValue( getRealPosition( position ));
    }

} // namespace Engine
