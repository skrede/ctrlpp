include(FetchContent)

FetchContent_Declare(tinympc
    GIT_REPOSITORY https://github.com/TinyMPC/TinyMPC.git
    GIT_TAG 023f36bd5b27267e1a96ed47a774d120b5707f13
    EXCLUDE_FROM_ALL)

FetchContent_MakeAvailable(tinympc)
