
def plotHistogram():
    ax.clear()
    ax.bar(R**3, N, width= binWidths, log=False, edgecolor="black", linewidth=1.5, color = colors, alpha=0.8)
    ax.set_xticks(R**3, labels=labels)
    ax.set_xlabel(fr"$(D [\mu m])^{{{3}}}$")
    ax.set_ylim([0,nDroplets*1.1])
    ax.vlines(R_eq**3, 0,1.2*np.max(N), colors='k', linestyles='--')
    ax.set_title(fr"$T = {T*1e3:.2f} ms \quad \alpha = {{{alpha:.2f}}}$")
    plt.draw()
    plt.pause(0.1)
    input("hit enter for next fig.")


def printLog(T, R,N):
    print(f"\tT = {1e3*T:.{3}f} ms")
    print(f"alpha = {calcAlpha(R,N):.{3}f}")
    printArrayWithNDecimals("N", N, 0)

def printArrayWithNDecimals(name, ar, nd):
    print(name, " = (", end="")
    for e in ar:
        print(f" {e:.{nd}f}", end=" ")
    print(f")\t\tsum = {np.sum(ar):.{nd}f}\n")
