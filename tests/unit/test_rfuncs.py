# import matplotlib.pyplot as plt
import numpy as np

from ciceroscm.carbon_cycle import rfuns


def test_rs_and_rb_functions_old():
    assert rfuns.rs_function(0) == 1.0
    assert rfuns.rs_function(3.0) == 0.6884435390654896
    assert rfuns.rb_function(0) == -3.599999999991197e-06


def test_rb_funcs():
    fld = rfuns.rb_function(np.arange(0, 10000))
    rbcoef = np.random.rand(3)
    rb_time = [
        5 * np.random.rand(1)[0],
        5 + 20 * np.random.rand(1)[0],
        25 + 50 * np.random.rand(1)[0],
    ]
    fld2 = rfuns.rb_function2(np.arange(0, 10000), rb_coef=rbcoef, rb_tim=rb_time)
    assert np.allclose(np.sum(fld), 24, rtol=5e-1)
    assert np.allclose(np.sum(fld2), 24, rtol=5.0e-1)
    assert np.allclose(fld[0], 0, atol=1e-5)
    assert np.allclose(fld2[0], 0, atol=1e-5)


def test_rs_funcs():
    fld = rfuns.rs_function_array(np.arange(0, 10000))
    rs_const = np.random.rand(1)[0] * 0.1
    rscoef = [
        rs_const,
        np.random.rand(1)[0],
        np.random.rand(1)[0],
        np.random.rand(1)[0],
    ]
    rs_time = [
        5 * np.random.rand(1)[0],
        5 + 20 * np.random.rand(1)[0],
        25 + 50 * np.random.rand(1)[0],
    ]
    fld2 = rfuns.rs_function2(np.arange(0, 10000), rs_coef=rscoef, rs_tim=rs_time)
    """
    fld3 = rfuns.rs_function2(np.arange(0,10000), rs_coef= [0.022936, 0.24278, 0.13963, 0.089318, 0.03782, 0.035549], rs_tim=[1.2679, 5.2528, 18.601, 68.736, 232.3])

    print(rs_time)
    print(rscoef)
    print(fld[:5])
    print(fld2[:5])
    print(fld3[:5])
    plt.plot(fld, label="orig")
    plt.plot(fld2, label = "upd")
    plt.plot(fld3, label = "long-term-orig")
    plt.legend()
    plt.show()
    print(np.sum([0.022936, 0.24278, 0.13963, 0.089318, 0.03782, 0.035549]))
    #assert np.allclose(np.sum(fld), 24, rtol=5e-2)
    #assert np.allclose(np.sum(fld2), 24, rtol= 5.e-2)
    #assert np.allclose(np.sum(fld3), np.sum(fld))
    """
    print(fld2[-1])
    print(rs_const / np.sum(rscoef))
    assert np.allclose(fld[0], 1, atol=1e-5)
    assert np.allclose(fld2[0], 1, atol=1e-5)
    assert np.allclose(fld2[-1], rs_const / np.sum(rscoef), atol=1e-2)
