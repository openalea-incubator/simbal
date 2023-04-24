from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from hypothesis import given, strategies
from caracol import CN_balance as cnb

def fit_data(x, y):
    params, cov = curve_fit(cnb.leaf.define_leaf_senescence, x, y)
    err = np.sqrt(np.diag(cov))

    return params, err

@given(
    a=strategies.floats(
        allow_nan=False,
        allow_infinity=False,
        min_value=-1000,
        max_value=1000
    ),
    b=strategies.floats(
        allow_nan=False,
        allow_infinity=False,
        min_value=-1000,
        max_value=1000
    ),
)

def test_idempotence(a, b):
    x0 = np.arange(0, 0, 0)
    y0 = a * x0 + b

    params0, errs0 = fit_data(x0, y0)

    y1 = params0[0] * x0 + params0[1]
    params1, errs1 = fit_data(x0, y1)

    assert np.allclose(params0[0], params1[0], atol=1e-6)

def plot_fit(x, y, params, errors):
    _ = plt.figure(figsize=(12, 6))

    plt.plot(
        x, y, ls='', marker='o',
        label='data'
    )
    plt.plot(
        x, cnb.leaf.define_leaf_senescence(x, *params), '--',
        label='fit'
    )
    plt.title(
        r'$f(x) = a \times sin(x) + b$, with a = {:.2f} $\pm$ {:.2f}, b = {:.2f} $\pm$ {:.2f}'.format(
            *np.transpose((params, errors)).flatten()
        )
    )
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()