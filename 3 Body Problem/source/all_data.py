import pandas as pd
from matplotlib import pyplot as plt

df = [pd.DataFrame() for i in range(15175)]

paths = []

n_phi = 0

def read_all():
    n_phi = 0
    phi = 0.400000000001
    while phi < 1.0:
        paths.append("./rsl/traj" + str(phi)[0:6] + ".txt")
        phi += 0.0001
        n_phi += 1

    print(n_phi)

    for i in range(n_phi):
        df[i] = pd.read_csv(paths[i], sep = " ", header = None)
        df[i].columns = ['Angle', 'Distance', 'Proximity', 'Time', 'Coll']
        if i%100 == 0: print(i)
    print('Read')
    return n_phi

def plot_1():
    df1 = pd.read_csv("./rsl/traj0.2052.txt", sep = " ", header = None)
    df1.columns = ['Angle', 'Distance', 'Proximty', 'Time', 'Coll']

    df1.plot(x = 'Time', y = 'Distance')
    plt.ylabel('Distance')

    plt.show()

def plot_all_distance():
    for i in range(1,n_phi):
        if df[i].Distance.min() < 0.01 and df[i].Distance.min() < df[i-1].Distance.min() and df[i].Distance.min() < df[i+1].Distance.min() or df[i].Distance.min() < 0.001:
            fig, ax = plt.subplots()
            ax.grid(True, color = 'gray', linestyle = 'dashed')
            #ax.title('Distance versus time')
            ax.set_xlabel('Time[s]')
            ax.set_ylabel('Distance')
            ax.plot(df[i].Time, df[i].Distance)
            fig.savefig('./rsl1/d'+str(i)+'.png')
            plt.close(fig)

def plot_all_proximity():
    for i in range(1,n_phi):
        if df[i].Distance.min() < 0.01 and df[i].Distance.min() < df[i-1].Distance.min() and df[i].Distance.min() < df[i+1].Distance.min() or df[i].Distance.min() < 0.001:
            fig, ax = plt.subplots()
            ax.grid(True, color = 'gray', linestyle = 'dashed')
            #ax.title('Distance versus time')
            ax.set_xlabel('Time[s]')
            ax.set_ylabel('Distance')
            ax.plot(df[i].Time, df[i].Proximity, color = 'red')
            fig.savefig('./rsl1/p'+str(i)+'.png')
            plt.close(fig)

n_phi = read_all()
plot_all_distance()
plot_all_proximity()
