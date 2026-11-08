//! Exercise golden comparison and blessing at the real workspace golden path.

#![forbid(unsafe_code)]

#[test]
fn real_workspace_golden_self_check() {
    iqsim_testkit::assert_golden("testkit/self_check", b"iqsim golden self-check\n");
}
