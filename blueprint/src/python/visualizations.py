import numpy as np
import matplotlib.pyplot as plt

def van_der_corput( k, alpha ):
    return max( alpha + (1 - k*alpha)/(2**k - 2), (1 - 2**(2-k))*alpha - (1-k*alpha)/(2**k - 2))

def optimized_van_der_corput( alpha ):
    bound = 1
    for k in range(2, 30):
        bound = min(bound, van_der_corput(k,alpha))
    return bound

def conjectured_bound(alpha):
    if alpha <= 1:
        return alpha/2
    else:
        return alpha-1

def van_der_corput_plot():
    alpha_range = np.linspace(0, 1, 500)
    plt.figure(figsize=(10, 6))

    plt.plot(alpha_range, [van_der_corput(2,alpha) for alpha in alpha_range], label=f'k=2 van der Corput')
    plt.plot(alpha_range, [van_der_corput(3,alpha) for alpha in alpha_range], label=f'k=3 van der Corput')
    plt.plot(alpha_range, [van_der_corput(4,alpha) for alpha in alpha_range], label=f'k=4 van der Corput')
    plt.plot(alpha_range, [van_der_corput(5,alpha) for alpha in alpha_range], label=f'k=5 van der Corput')
    plt.plot(alpha_range, [optimized_van_der_corput(alpha) for alpha in alpha_range], label=f'Optimized van der Corput')
    plt.plot(alpha_range, [alpha/2 for alpha in alpha_range], label=f'Conjectured value')
    plt.xlabel('alpha')
    plt.ylabel('beta(alpha)') 
    plt.title('Classical bounds on beta(alpha)')
    plt.legend()
    plt.grid(True)
    plt.ylim(0, 1)
    plt.show()

def van_der_corput_plot2():
    alpha_range = np.linspace(0, 1, 500)
    alpha_full_range = np.linspace(0, 2, 1000)

    plt.figure(figsize=(10, 6))

    plt.plot(alpha_full_range, [optimized_van_der_corput(alpha) for alpha in alpha_full_range], label=f'Optimized van der Corput')
    plt.plot(alpha_range, [alpha for alpha in alpha_range], label=f'Trivial upper bound')
    plt.plot(alpha_full_range, [conjectured_bound(alpha) for alpha in alpha_full_range], label=f'Conjectured value')
    plt.xlabel('alpha')
    plt.ylabel('beta(alpha)') 
    plt.title('Classical bounds on beta(alpha)')
    plt.legend()
    plt.grid(True)
    plt.ylim(0, 1)
    plt.show()

van_der_corput_plot2()