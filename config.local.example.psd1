@{
    # Copy to config.local.psd1 (ignored by Git); see docs/setup.md.
    Sbcl = 'sbcl' # Supported Windows x86-64 installation.
    Python = 'python' # Python 3.11+; standard library only.
    CoaltonRoot = '' # Absolute installed checkout matching toolchain.lock.json.
    QuicklispRoot = '' # Absolute isolated Quicklisp root matching the lock.
    CacheRoot = '' # Empty uses build/cache. Use a separate cache per checkout; a new directory for a fresh build.
}