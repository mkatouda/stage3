set(GENERATE_DOCS TRUE)

find_package(Sphinx QUIET)
find_package(LATEX QUIET)

if(NOT SPHINX_FOUND)
    message(WARNING "Cannot find Sphinx. Cannot generate documentation")
    set(GENERATE_DOCS FALSE)
endif()

if(GENERATE_DOCS)

    if(NOT DEFINED SPHINX_THEME)
        set(SPHINX_THEME default)
    endif()

    if(NOT DEFINED SPHINX_THEME_DIR)
        set(SPHINX_THEME_DIR)
    endif()

    # configured documentation tools and intermediate build results
    set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")

    # Sphinx cache with pickled ReST documents
    set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

    # HTML output directory
    set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")

    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in"
        "${BINARY_BUILD_DIR}/conf.py"
        @ONLY)

    add_custom_target(docs_html ALL
        ${SPHINX_EXECUTABLE}
            -q -b html
            -c "${BINARY_BUILD_DIR}"
            -d "${SPHINX_CACHE_DIR}"
            "${CMAKE_CURRENT_SOURCE_DIR}"
            "${SPHINX_HTML_DIR}"
        COMMENT "Building HTML documentation with Sphinx")

    list(APPEND CLEAN_FILES "${SPHINX_CACHE_DIR}" "${SPHINX_HTML_DIR}")

    add_dependencies(docs_html docs_sources)

    if(PDFLATEX_COMPILER)
        # PDF output directory
        set(SPHINX_PDF_DIR "${CMAKE_CURRENT_BINARY_DIR}/pdf")

        add_custom_target(docs_pdf ALL
            ${SPHINX_EXECUTABLE}
                -q -b latex
                -c "${BINARY_BUILD_DIR}"
                -d "${SPHINX_CACHE_DIR}"
                "${CMAKE_CURRENT_SOURCE_DIR}"
                "${SPHINX_PDF_DIR}"
            COMMENT "Building PDF documentation with Sphinx")

        add_custom_command(TARGET docs_pdf POST_BUILD
                           COMMAND ${CMAKE_MAKE_PROGRAM} LATEXOPTS=-interaction=batchmode
                           WORKING_DIRECTORY "${SPHINX_PDF_DIR}")

        list(APPEND CLEAN_FILES "${SPHINX_PDF_DIR}")

        add_dependencies(docs_pdf docs_sources)

    else()
        message(WARNING "Cannot find LaTeX for build PDF documentation")
    endif()

    # Add output directories to clean target.
    set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEAN_FILES}")


endif(GENERATE_DOCS)