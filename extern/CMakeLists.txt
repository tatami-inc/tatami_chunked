cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

include(FetchContent)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/tatami-inc/tatami
  GIT_TAG aa8a89545999826741d3225f0b2fbdf62b6bf52c
)

FetchContent_MakeAvailable(tatami)
