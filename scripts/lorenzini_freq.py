"""Frequency tuning of a Lorenzini ampulla: the dielectric (capacitive) part.

DC conduction gives the sensitivity but not the well-known band-pass tuning (~0.1 to 10 Hz)
of the ampullae. Adding the dielectric, kappa = sigma + j omega eps, the response becomes
frequency dependent. In the equivalent circuit the gel canal is a series access resistance
R_a (which our field model extracts), the sensory epithelium is a leaky capacitor R_m // C_m,
and the surrounding skin / canal wall shunts the very low frequencies. The product of a
high-pass (skin/wall) and a low-pass (membrane) is a band-pass, exactly the measured tuning.
We compute the transfer H(omega) = V_membrane / V_external and locate the pass band.
"""
import os
import numpy as np

# Figure output directory. Defaults to the article folder, overridable so the repository
# runs anywhere (a hard-coded absolute path is not reproducible for a third party).
FIGDIR = os.environ.get("LORENZINI_FIGDIR", os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures"))
os.makedirs(FIGDIR, exist_ok=True)


# access (canal) resistance, from the field model: R_a = L /(sigma_gel * A_canal)
L_CAN, R_CAN, SIG_GEL = 0.05, 0.5e-3, 4.0
R_A = L_CAN / (SIG_GEL * np.pi * R_CAN**2)

# sensory epithelium (leaky capacitor, low-pass) and skin / canal-wall (high-pass, adaptation)
R_M, C_M = 6.6e6, 3.0e-9        # membrane RC -> low-pass corner ~8 Hz
R_W, C_W = 2.0e6, 0.4e-6        # skin / canal-wall RC -> high-pass corner ~0.2 Hz

f1 = 1 / (2*np.pi*R_W*C_W)      # high-pass corner (Hz)
f2 = 1 / (2*np.pi*R_M*C_M)      # low-pass corner (Hz)


def H(f):
    """Band-pass = (skin/wall high-pass) x (membrane low-pass). The canal access R_A sets the
    absolute sensitivity (from the DC field model); the two RC corners set the pass band."""
    w = 2*np.pi*f
    hp = (1j*w*R_W*C_W) / (1 + 1j*w*R_W*C_W)
    lp = 1.0 / (1 + 1j*w*R_M*C_M)
    return hp * lp


def figure():
    global C_M
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    f = np.geomspace(1e-2, 1e3, 400)
    Hf = np.abs(H(f)); Hn = Hf / Hf.max()
    fpk = f[np.argmax(Hf)]
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    ax[0].semilogx(f, 20*np.log10(Hn), color="#185FA5", lw=2)
    ax[0].axvspan(f1, f2, color="#cfe6cf", alpha=.6, label=f"pass band ~{f1:.1f}-{f2:.0f} Hz")
    ax[0].axvline(fpk, color="#A32D2D", ls="--", lw=1, label=f"peak {fpk:.1f} Hz")
    ax[0].axhline(-3, color="0.6", lw=.8, ls=":")
    ax[0].set_xlabel("frequency (Hz)"); ax[0].set_ylabel("normalised response (dB)")
    ax[0].set_ylim(-40, 2); ax[0].legend(fontsize=9, loc="lower center"); ax[0].grid(alpha=.3, which="both")
    ax[0].set_title("(a) Ampulla band-pass tuning from the dielectric\n(membrane capacitance + skin/wall shunt)")
    tau_w, tau_m = R_W*C_W, R_M*C_M
    ax[0].text(0.03, 0.97,
               r"$H(\omega)=\dfrac{j\omega\tau_w}{1+j\omega\tau_w}\cdot\dfrac{1}{1+j\omega\tau_m}$"
               "\n"
               r"$f_1=\dfrac{1}{2\pi\tau_w},\ f_2=\dfrac{1}{2\pi\tau_m},\ f_0=\dfrac{1}{2\pi\sqrt{\tau_w\tau_m}}$"
               "\n"
               fr"$f_1\!\approx\!{f1:.1f},\ f_0\!\approx\!{fpk:.1f},\ f_2\!\approx\!{f2:.0f}$ Hz",
               transform=ax[0].transAxes, va="top", ha="left", fontsize=8.5,
               bbox=dict(boxstyle="round", fc="#fbf7ee", ec="0.6"))
    # corner sensitivity to membrane capacitance
    for cm, c in [(1e-9, "#888780"), (3e-9, "#185FA5"), (1e-8, "#A32D2D")]:
        C_M = cm
        Hf2 = np.abs(H(f))
        ax[1].semilogx(f, Hf2/np.abs(H(f)).max() if False else Hf2/Hf2.max(), color=c, lw=1.8, label=f"C_m={cm*1e9:.0f} nF")
    C_M = 3e-9
    ax[1].set_xlabel("frequency (Hz)"); ax[1].set_ylabel("normalised response")
    ax[1].set_title("(b) The membrane capacitance sets the upper corner\n(tunable low-pass edge)")
    ax[1].legend(fontsize=9); ax[1].grid(alpha=.3, which="both")
    fig.suptitle("Dielectric (frequency) response of a Lorenzini ampulla: a certified band-pass electroreceptor", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "Fig_lorenzini_freq.png"), dpi=125)
    print(f"R_a (canal access) = {R_A/1e3:.1f} kOhm")
    print(f"high-pass corner f1 = {f1:.2f} Hz ; low-pass corner f2 = {f2:.1f} Hz ; peak ~ {fpk:.1f} Hz")
    print("-> saved Fig_lorenzini_freq.png")


if __name__ == "__main__":
    figure()
