# repulsive
# To avoid copying things to GPU memory,
# ideally allocate everything in torch on the GPU
# and avoid non-torch function calls
import torch
from torchquad import VEGAS, set_up_backend, MonteCarlo
from itertools import permutations
import datetime

# flake8: noqa 
# noqa: E501

torch._dynamo.config.cache_size_limit = 64

set_up_backend("torch", data_type="float32")
base_num_ints = 10000000
num_ints = base_num_ints

display_plots = False

def wf_single(x):
    return 1.0-torch.exp(-(x[:, 0]-x[:, 1])**2)

factor = 5.0


def wf(x):
    # first half of x is x0, second half is x1 (two different atoms)
    ret = 1.0
    for i in range(0, x.shape[1] // 2):
        ret = ret * (1-torch.exp(-(torch.fmod(x[:, i]-x[:, i+x.shape[1] // 2] + factor*torch.pi, 2*factor*torch.pi) - factor*torch.pi)**2))
    return ret


def wf_sym(x):
    def perm_parity(lst):
        '''\
        Given a permutation of the digits 0..N in order as a list,
        returns its parity (or sign): +1 for even parity; -1 for odd.
        '''
        parity = 1
        for i in range(0, len(lst)-1):
            if lst[i] != i:
                parity *= -1
                mn = min(range(i, len(lst)), key=lst.__getitem__)
                lst[i], lst[mn] = lst[mn], lst[i]
        return parity

    dd = x.shape[1] // 2
    dindicies = range(dd)
    pdindices1 = permutations(dindicies)
    pdindices2 = permutations(dindicies)
    ret = 0.0
    x0 = x[:, :dd]
    x1 = x[:, dd:]
    for a in pdindices1:
        p1 = perm_parity(list(a))
        for b in pdindices2:
            ret = ret + p1 * perm_parity(list(b))*wf(torch.cat((x0[:, a], x1[:, b]), dim=1))
    return ret


def rho(x):
    return wf_sym(x)**2


# Declare an integrator;
# here we use the simple, stochastic Monte Carlo integration method
mc = MonteCarlo()

print("integrator is ", mc)

def rho_od(a1s, b1s, a1, b1):
    a1s = torch.tensor([[a1s]])
    b1s = torch.tensor([[b1s]])
    a1 = torch.tensor([[a1]])
    b1 = torch.tensor([[b1]])
    with torch.no_grad():
        integral_value = mc.integrate(
                        torch.compile(lambda x: wf_sym(torch.cat([a1s.expand((x.shape[0],1)), x[:, :dd], b1s.expand((x.shape[0],1)), x[:, dd:]],dim=1)) *
                        wf_sym(torch.cat([a1.expand((x.shape[0],1)), x[:, :dd], b1.expand((x.shape[0],1)), x[:, dd:]],dim=1))),
                        dim=dd,
                        N=num_ints,
                        integration_domain= dd * [[-factor*torch.pi, factor*torch.pi]],
                    )
    return integral_value

dd = 4

def rhotest(x):
    return rho_od(x, x, 0.0, 0.0)/torch.sqrt(rho_od(x, x, x, x)*rho_od(0.0, 0.0, 0.0, 0.0))


print(datetime.datetime.now(), "Integral value: ", rhotest(1.2))
print(datetime.datetime.now())

import numpy as np
from matplotlib import pyplot as plt

plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True


do_plots = False
if do_plots:
    x = np.linspace(-factor*3.0, factor*3.0, 50)

    dd = 4
    y = []
    for i in x:
        y.append(rhotest(float(i)))
    plt.plot(x,y)

    plt.savefig("repulsive_"+str(dd)+".png")

    dd = 6
    y = []
    for i in x:
        y.append(rhotest(float(i)))
    plt.plot(x,y)

    plt.savefig("repulsive_"+str(dd)+".png")

    dd = 8
    y = []
    for i in x:
        y.append(rhotest(float(i)))
    plt.plot(x,y)

    plt.savefig("repulsive_"+str(dd)+".png")

    dd = 10
    y = []
    for i in x:
        y.append(rhotest(float(i)))
    plt.plot(x,y)

    # plt.savefig("repulsive_"+str(dd)+".png")

    # dd = 12
    # y = []
    # for i in x:
    #     y.append(rhotest(float(i)))
    # plt.plot(x,y)
    plt.title("Repulsive factor="+str(factor))
    plt.xlabel("x")
    plt.ylabel("correlation coefficient")
    plt.savefig("repulsive_all.png")

    if display_plots:
        plt.show()
    else:
        plt.clf()

do_calcs = True
if do_calcs:
    num_ints = 5 * base_num_ints
    dds = [4,6,8,10,12]
    y = []
    for x in dds:
        dd = x
        y.append(rhotest(10).cpu())
        print(datetime.datetime.now(), dd, y[-1])
        plt.plot([i + 2 for i in dds[:len(y)]], y)
        plt.title("Repulsive factor="+str(factor)+" at x=10")
        plt.xlabel("number of particles")
        plt.ylabel("correlation coefficient")
        plt.savefig("repulsive_dd.png")
    if display_plots:
        plt.show()
    else:
        plt.clf()

do_fact = True
if do_fact:
    dd = 4
    num_ints = 5 * base_num_ints
    dds = [1,2,3,4,5,6,7,8,9,10]
    y = []
    for x in dds:
        factor = x
        y.append(rhotest(2*factor).cpu())
        print(datetime.datetime.now(), x, y[-1])
        plt.plot(dds[:len(y)], y)
        plt.title("Repulsive particles="+str(dd +2 )+" at x=2*factor")
        plt.xlabel("factor of the x periode")
        plt.ylabel("correlation coefficient")
        plt.savefig("repulsive_factor.png")
    if display_plots:
        plt.show()
    else:
        plt.clf()

do_scale = True
if do_scale:
    num_ints = 5 * base_num_ints
    dds = [2,4,6,8,10]
    y = []
    for x in dds:
        dd = x
        factor = (x + 2) // 2
        y.append(rhotest(2*factor).cpu())
        print(datetime.datetime.now(), dd, y[-1])
        plt.plot([i+2 for i in dds[:len(y)]], y)
        plt.title("Repulsive particles at x=2*factor")
        plt.xlabel("number of particles with factor=1/2*number of particles")
        plt.ylabel("correlation coefficient")
        plt.savefig("repulsive_scale.png")
    if display_plots:
        plt.show()
    else:
        plt.clf()
