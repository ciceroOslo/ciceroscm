from ciceroscm import _utils


def test_check_numeric_pamset():
    required = {"Year": 2022, "number": 4.55, "test": 6.3e-4, "test2": 85}
    pamset = {
        "test": 7.8e-3,
        "test2": "string",
        "untouched": 44,
        "untouched_string": "Hello",
    }
    pamset2 = _utils.check_numeric_pamset(required, pamset)
    assert pamset2["Year"] == 2022
    assert pamset2["number"] == 4.55
    assert pamset2["test"] == 7.8e-3
    assert pamset2["test2"] == 85
    assert pamset2["untouched"] == 44
    assert pamset2["untouched_string"] == "Hello"


def test_cut_and_check_pamset():
    required = {"Year": 2022, "number": 4.55, "test": 6.3e-4, "test2": 85}
    pamset = {
        "test": 7.8e-3,
        "test2": "string",
        "untouched": 44,
        "untouched_string": "Hello",
    }
    used = {
        "test": 6,
        "test2": "string",
        "untouched": 44,
    }
    pamset2 = _utils.cut_and_check_pamset(
        required, pamset, used=used, cut_warnings=True
    )
    assert pamset2["Year"] == 2022
    assert pamset2["number"] == 4.55
    assert pamset2["test"] == 7.8e-3
    assert pamset2["test2"] == 85
    expected_keys = ["test", "test2", "untouched", "Year", "number"]
    for key in expected_keys:
        assert key in pamset2
    assert "untouched_string" not in pamset2


def test_update_pam_if_numeric():
    pamset_orig = {
        "beta_f": 0.287,
        "other_pam": 1,
        "nystart": 1750,
        "leave_it": "Seriously, leave it",
    }

    pamset_new = {"beta_f": 3, "weird": 43, "other_pam": "hello", "nystart": 1850}
    pamset_expect = {
        "beta_f": 3,
        "other_pam": 1,
        "nystart": 1750,
        "leave_it": "Seriously, leave it",
    }
    assert pamset_orig == _utils.update_pam_if_numeric(pamset_orig.copy(), None, None)
    assert pamset_expect == _utils.update_pam_if_numeric(
        pamset_orig.copy(), pamset_new, ["beta_f", "other_pam", "weird"]
    )
