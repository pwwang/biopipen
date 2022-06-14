patch \
    $(python -c 'from poetry.repositories import installed_repository as ir; print(ir.__file__, end="")') \
    tests/conda/poetry.pth
