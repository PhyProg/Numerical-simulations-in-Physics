import numpy as np
import matplotlib.pyplot as plt

from density_calculator import DensityOfStatesCalculator
from dispersion import quadratic_lattice_calculator, graphene_dispersion_calculator


def plot_quadratic_lattice(no_of_grid_points = 200, no_of_bins = 'auto'):
    disp = quadratic_lattice_calculator()
    grid = np.array([np.linspace(-np.pi, np.pi, no_of_grid_points), np.linspace(-np.pi, np.pi, no_of_grid_points)])

    calc = DensityOfStatesCalculator(disp, grid)
    
    fig, ax = calc.get_3d_plot()

    fig.set_figheight(9)
    fig.set_figwidth(15)

    ax.set_title("Dispersion relation for quadratic lattice")

    ax.set_zlabel(r"$\frac{\epsilon_k - \epsilon_0}{2t}$")
    ax.set_ylabel(r"$k_y a_y$")
    ax.set_xlabel(r"$k_x a_x$")

    fig.savefig("../plot/quadratic_3d.png")

    fig, ax = calc.get_2d_contourf_plot()

    fig.savefig("../plot/quadratic_2d.png")

    fig, ax = calc.get_histogram_plot(no_of_bins)
    
    fig.savefig("../plot/quadratic_hist.png")



def plot_honeycomb_lattice(no_of_grid_points = 200, no_of_bins = 'auto'):
    disp = graphene_dispersion_calculator()
    grid = np.array([np.linspace(-np.pi, np.pi, no_of_grid_points), np.linspace(-np.pi, np.pi, no_of_grid_points)])

    calc = DensityOfStatesCalculator(disp, grid)
    
    fig, ax = calc.get_3d_plot()

    fig.set_figheight(9)
    fig.set_figwidth(15)

    ax.set_title("Dispersion relation for honeycomb lattice (Graphene)")

    ax.set_zlabel(r"$\frac{\epsilon_k - \epsilon_0}{2t}$")
    ax.set_ylabel(r"$k_y a_y$")
    ax.set_xlabel(r"$k_x a_x$")

    fig.savefig("../plot/graphene_3d.png")

    fig, ax = calc.get_2d_contourf_plot()

    fig.savefig("../plot/graphene_2d.png")

    fig, ax = calc.get_histogram_plot(no_of_bins)
    
    fig.savefig("../plot/honeycomb_hist.png")


if __name__ == "__main__":
    plot_quadratic_lattice(no_of_grid_points = 6000, no_of_bins = 100000)
    plot_honeycomb_lattice(no_of_grid_points = 6000, no_of_bins = 100000)

