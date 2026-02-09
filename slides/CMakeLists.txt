find_program(MARP_EXECUTABLE marp REQUIRED)

set(SLIDES_MD ${CMAKE_CURRENT_SOURCE_DIR}/slides.md)
set(SLIDES_HTML ${CMAKE_BINARY_DIR}/index.html)

add_custom_command(
        OUTPUT ${SLIDES_HTML}
        COMMAND ${MARP_EXECUTABLE}
        ${SLIDES_MD}
        --html
        --output ${SLIDES_HTML}
        DEPENDS ${SLIDES_MD}
        COMMENT "Generating Marp HTML slides"
        VERBATIM
)

add_custom_target(slides ALL
        DEPENDS ${SLIDES_HTML}
)
