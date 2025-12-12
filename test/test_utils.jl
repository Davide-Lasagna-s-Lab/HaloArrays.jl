import HaloArrays: swapregions,
                   LEFT, RIGHT, ALL, CENTER,
                   subarray_slices,
                   SEND, RECV

using Test

@testset "swapregions                            " begin
    # 1d
    topology    = (  10,)
    isperiodic  = (true,)
    @test swapregions(topology, isperiodic, true) == [(LEFT,), (RIGHT,)]

    topology    = (   10,)
    isperiodic  = (false,)
    @test swapregions(topology, isperiodic, true) == [(LEFT,), (RIGHT,)]

    topology    = (   1,)
    isperiodic  = (true,)
    @test swapregions(topology, isperiodic, true) == [(LEFT,), (RIGHT,)]

    topology    = (    1,)
    isperiodic  = (false,)
    @test swapregions(topology, isperiodic, true) == []

    # 2d
    topology    = (   1,     1)
    isperiodic  = (true, false)
    @test swapregions(topology, isperiodic, true) == [(LEFT,   CENTER),
                                                      (RIGHT,  CENTER)]

    topology    = (   10,   10)
    isperiodic  = (false, false)
    @test swapregions(topology, isperiodic, true) == [(LEFT,   CENTER),
                                                      (RIGHT,  CENTER),
                                                      (CENTER, LEFT),
                                                      (CENTER, RIGHT)]

    topology    = (  2,     2)
    isperiodic  = (true, true)
    @test swapregions(topology, isperiodic, true) == [(LEFT,   CENTER),
                                                      (RIGHT,  CENTER),
                                                      (CENTER, LEFT),
                                                      (CENTER, RIGHT)]

    topology    = (  1,     1)
    isperiodic  = (true, true)
    @test swapregions(topology, isperiodic, true) == [(LEFT,   CENTER),
                                                      (RIGHT,  CENTER),
                                                      (CENTER, LEFT),
                                                      (CENTER, RIGHT)]

    topology    = (  1,      1)
    isperiodic  = (false, true)
    @test swapregions(topology, isperiodic, true) == [(CENTER, LEFT),
                                                      (CENTER, RIGHT)]

    # 3d
    topology    = (    1,     1,     1)
    isperiodic  = (false, false, false)
    @test swapregions(topology, isperiodic, true) == []

    topology    = (  1,      10,     1)
    isperiodic  = (false, false, false)
    @test swapregions(topology, isperiodic, true) == [(CENTER,  LEFT, CENTER),
                                                      (CENTER, RIGHT, CENTER)]

    topology    = (  1,      10,    10)
    isperiodic  = (false, false, false)
    @test swapregions(topology, isperiodic, true) == [(CENTER,   LEFT, CENTER),
                                                      (CENTER,  RIGHT, CENTER),
                                                      (CENTER, CENTER, LEFT),
                                                      (CENTER, CENTER, RIGHT)]

    topology    = (   1,    1,     10)
    isperiodic  = (true, false, false)
    @test swapregions(topology, isperiodic, true) == [(  LEFT, CENTER, CENTER),
                                                      ( RIGHT, CENTER, CENTER),
                                                      (CENTER, CENTER, LEFT),
                                                      (CENTER, CENTER, RIGHT)]

    topology    = (   1,    10,    1)
    isperiodic  = (true, false, true)
    @test swapregions(topology, isperiodic, true) == [(LEFT,   CENTER, CENTER),
                                                      (RIGHT,  CENTER, CENTER),
                                                      (CENTER, LEFT,   CENTER),
                                                      (CENTER, RIGHT,  CENTER),
                                                      (CENTER, CENTER, LEFT),
                                                      (CENTER, CENTER, RIGHT)]

    topology    = (   1,    10,    1)
    isperiodic  = (true, false, true)
    @test swapregions(topology, isperiodic, false) == [(LEFT,  ALL,   ALL),
                                                       (RIGHT, ALL,   ALL),
                                                       (ALL,   LEFT,  ALL),
                                                       (ALL,   RIGHT, ALL),
                                                       (ALL,   ALL,   LEFT),
                                                       (ALL,   ALL,   RIGHT)]
end

@testset "subarray slices                        " begin
    # 1d

    # 1 2 3 4 5 6 7 8 9
    # h h 1 2 3 4 5 h h
    @test subarray_slices((5,), (2,), (LEFT,),   SEND) == (3:4,)
    @test subarray_slices((5,), (2,), (RIGHT,),  SEND) == (6:7,)
    @test subarray_slices((5,), (2,), (CENTER,), SEND) == (3:7,)
    @test subarray_slices((5,), (2,), (ALL,),    SEND) == (1:9,)
    @test subarray_slices((5,), (2,), (LEFT,),   RECV) == (1:2,)
    @test subarray_slices((5,), (2,), (RIGHT,),  RECV) == (8:9,)
    @test subarray_slices((5,), (2,), (CENTER,), RECV) == (3:7,)
    @test subarray_slices((5,), (2,), (ALL,),    RECV) == (1:9,)

    # 1 2 3 4 5 6 7 8 9 10 11
    # h h h 1 2 3 4 5 h  h  h
    @test subarray_slices((5,), (3,), (LEFT,),   SEND) == (4:6,)
    @test subarray_slices((5,), (3,), (RIGHT,),  SEND) == (6:8,)
    @test subarray_slices((5,), (3,), (CENTER,), SEND) == (4:8,)
    @test subarray_slices((5,), (3,), (ALL,),    SEND) == (1:11,)
    @test subarray_slices((5,), (3,), (LEFT,),   RECV) == (1:3,)
    @test subarray_slices((5,), (3,), (RIGHT,),  RECV) == (9:11,)
    @test subarray_slices((5,), (3,), (CENTER,), RECV) == (4:8,)
    @test subarray_slices((5,), (3,), (ALL,),    RECV) == (1:11,)

    # 2d
    @test subarray_slices((5, 4), (2, 2), (LEFT,   ALL),    SEND) == (3:4, 1:8)
    @test subarray_slices((5, 4), (2, 2), (LEFT,   CENTER), SEND) == (3:4, 3:6)
    @test subarray_slices((5, 4), (2, 2), (CENTER, LEFT),   SEND) == (3:7, 3:4)
    @test subarray_slices((5, 4), (2, 2), (ALL,    RIGHT),  SEND) == (1:9, 5:6)
    @test subarray_slices((5, 4), (2, 2), (LEFT,   ALL),    RECV) == (1:2, 1:8)
    @test subarray_slices((5, 4), (2, 2), (LEFT,   CENTER), RECV) == (1:2, 3:6)
    @test subarray_slices((5, 4), (2, 2), (CENTER, LEFT),   RECV) == (3:7, 1:2)
    @test subarray_slices((5, 4), (2, 2), (ALL,    RIGHT),  RECV) == (1:9, 7:8)

    # 3d
    @test subarray_slices((5, 4, 3), (1, 2, 3), (LEFT,   ALL,    ALL),    SEND) == (2:2, 1:8, 1:9)
    @test subarray_slices((5, 4, 3), (1, 2, 3), (RIGHT,  ALL,    ALL),    SEND) == (6:6, 1:8, 1:9)
    @test subarray_slices((5, 4, 3), (1, 2, 3), (LEFT,   CENTER, CENTER), SEND) == (2:2, 3:6, 4:6)
    @test subarray_slices((5, 4, 3), (1, 2, 3), (CENTER, LEFT,   ALL),    SEND) == (2:6, 3:4, 1:9)
    @test subarray_slices((5, 4, 3), (1, 2, 3), (CENTER, RIGHT,  ALL),    SEND) == (2:6, 5:6, 1:9)
    @test subarray_slices((5, 4, 3), (1, 2, 3), (CENTER, LEFT,   CENTER), SEND) == (2:6, 3:4, 4:6)
    @test subarray_slices((5, 4, 3), (1, 2, 3), (CENTER, RIGHT,  CENTER), SEND) == (2:6, 5:6, 4:6)
    @test subarray_slices((5, 4, 3), (1, 2, 3), (CENTER, CENTER, RIGHT),  SEND) == (2:6, 3:6, 4:6)
    @test subarray_slices((5, 4, 3), (1, 2, 3), (LEFT,   ALL,    ALL),    RECV) == (1:1, 1:8, 1:9)
    @test subarray_slices((5, 4, 3), (1, 2, 3), (RIGHT,  ALL,    ALL),    RECV) == (7:7, 1:8, 1:9)
    @test subarray_slices((5, 4, 3), (1, 2, 3), (LEFT,   CENTER, CENTER), RECV) == (1:1, 3:6, 4:6)
    @test subarray_slices((5, 4, 3), (1, 2, 3), (CENTER, LEFT,   ALL),    RECV) == (2:6, 1:2, 1:9)
    @test subarray_slices((5, 4, 3), (1, 2, 3), (CENTER, RIGHT,  ALL),    RECV) == (2:6, 7:8, 1:9)
    @test subarray_slices((5, 4, 3), (1, 2, 3), (CENTER, LEFT,   CENTER), RECV) == (2:6, 1:2, 4:6)
    @test subarray_slices((5, 4, 3), (1, 2, 3), (CENTER, RIGHT,  CENTER), RECV) == (2:6, 7:8, 4:6)
    @test subarray_slices((5, 4, 3), (1, 2, 3), (CENTER, CENTER, RIGHT),  RECV) == (2:6, 3:6, 7:9)
    # these should be enough. the algorithm works for any number of dimensions
end
