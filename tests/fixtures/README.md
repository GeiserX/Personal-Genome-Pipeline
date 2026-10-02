# Test fixture

`VERSION` names the GitHub release that holds the e2e test data: a small slice
of the public GIAB HG002 genome, built by `scripts/ci/build-fixture.sh`. The
data itself is not in git.

To change the data, edit the build script and bump `VERSION` in the same
commit; a push builds and publishes the new release. Contents, regions and
the full procedure: [docs/testing.md](../../docs/testing.md).
