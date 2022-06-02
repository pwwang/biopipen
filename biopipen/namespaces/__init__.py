import warnings

warnings.simplefilter('always', DeprecationWarning)
warnings.warn(
    "`biopipen.namespaces` is deprecated. Please use `biopipen.ns` instead.",
    DeprecationWarning,
)
