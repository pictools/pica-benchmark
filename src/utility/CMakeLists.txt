include_directories(
    ../
    ../../pica/include
    ../../thirdParty)

add_library(benchmark-utility
    FieldGenerator.h
    Output.h
    Parameters.h
    Parser.h
    ParticleGenerator.h
    Random.h
    Timer.h

    Output.cpp
    Parser.cpp
    Random.cpp
    Timer.cpp)
