import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

def set_styles():
    plt.style.use(['seaborn-white', 'seaborn-paper'])
#     matplotlib.rc("font", family="Times New Roman")

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42