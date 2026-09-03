import numpy as np
import matplotlib.pyplot as plt
import os

# import sys
# sys.path.append("C:\\Program Files\\Lumerical\\v241\\api\\python")
import lumapi

# 1. Setup Parameters
start_wl_nm = 1530
stop_wl_nm = 1570
c = 299792458
start_freq = c / (stop_wl_nm * 1e-9)
stop_freq = c / (start_wl_nm * 1e-9)
dat_file_path = os.path.abspath("psgc_lumerical.dat")

# 2. Launch Session
print("Launching Lumerical INTERCONNECT...")
intc = lumapi.INTERCONNECT(hide=False)
intc.switchtodesign()
intc.deleteall()

# 3. Create PSGC
intc.addelement("Optical N Port S-Parameter")
intc.set("name", "PSGC")
intc.set("load from file", True)
intc.set("s parameters filename", dat_file_path)

# 4. Create a SINGLE ONA
intc.addelement("Optical Network Analyzer")
intc.set("name", "ONA")
intc.set("input parameter", 2)
intc.set("start frequency", start_freq)
intc.set("stop frequency", stop_freq)
intc.set("number of points", 1000)
intc.set("plot kind", "wavelength")

# =========================================================================
# PHASE 1: FORWARD INJECTION (Fiber -> Chip)
# =========================================================================
intc.setnamed("ONA", "number of input ports", 2)

# Wire up the forward connections
intc.connect("ONA", "output", "PSGC", "Fiber")
intc.connect("PSGC", "WG_X", "ONA", "input 1")
intc.connect("PSGC", "WG_Y", "ONA", "input 2")

# --- Test 1: Fiber TE (X-Pol) ---
print("Running Fiber TE Injection...")
intc.setnamed("ONA", "orthogonal identifier", 1)
intc.run()

ds_F_TE_1 = intc.getresult("ONA", "input 1/mode 1/gain")
ds_F_TE_2 = intc.getresult("ONA", "input 2/mode 1/gain")
wavelengths_nm = ds_F_TE_1["wavelength"].flatten() * 1e9  # Save shared X-axis

F_TE_to_WG_X = ds_F_TE_1["TE gain (dB)"].flatten()
F_TE_to_WG_Y = ds_F_TE_2["TE gain (dB)"].flatten()

# --- Test 2: Fiber TM (Y-Pol) ---
print("Running Fiber TM Injection...")
intc.switchtodesign()
intc.setnamed("ONA", "orthogonal identifier", 2)
intc.run()

ds_F_TM_1 = intc.getresult("ONA", "input 1/mode 1/gain")
ds_F_TM_2 = intc.getresult("ONA", "input 2/mode 1/gain")

F_TM_to_WG_X = ds_F_TM_1["TE gain (dB)"].flatten()
F_TM_to_WG_Y = ds_F_TM_2["TE gain (dB)"].flatten()

# =========================================================================
# PHASE 2: REVERSE INJECTION (Chip -> Fiber)
# =========================================================================
intc.switchtodesign()

# Tear down forward connections
intc.disconnect("ONA", "output", "PSGC", "Fiber")
intc.disconnect("PSGC", "WG_X", "ONA", "input 1")
intc.disconnect("PSGC", "WG_Y", "ONA", "input 2")

# Reconfigure ONA to only measure the Fiber output port
intc.setnamed("ONA", "number of input ports", 1)
intc.connect("PSGC", "Fiber", "ONA", "input 1") # The Fiber output stays plugged into Input 1

# --- Test 3: Waveguide X Injection ---
print("Running WG_X Injection...")
intc.connect("ONA", "output", "PSGC", "WG_X")
intc.run()

ds_X_m1 = intc.getresult("ONA", "input 1/mode 1/gain") # Captures Fiber TE
ds_X_m2 = intc.getresult("ONA", "input 1/mode 2/gain") # Captures Fiber TM

WG_X_to_F_TE = ds_X_m1["TE gain (dB)"].flatten()
WG_X_to_F_TM = ds_X_m2["TE gain (dB)"].flatten()

# --- Test 4: Waveguide Y Injection ---
print("Running WG_Y Injection...")
intc.switchtodesign()
intc.disconnect("ONA", "output", "PSGC", "WG_X")
intc.connect("ONA", "output", "PSGC", "WG_Y")
intc.run()

ds_Y_m1 = intc.getresult("ONA", "input 1/mode 1/gain")
ds_Y_m2 = intc.getresult("ONA", "input 1/mode 2/gain")

WG_Y_to_F_TE = ds_Y_m1["TE gain (dB)"].flatten()
WG_Y_to_F_TM = ds_Y_m2["TE gain (dB)"].flatten()

print("Simulations complete! Generating plots...")

# =========================================================================
# MATPLOTLIB DASHBOARD
# =========================================================================
fig, axs = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('PSGC Bidirectional S-Parameter Verification', fontsize=16)

# Plot 1: Forward TE
axs[0, 0].plot(wavelengths_nm, F_TE_to_WG_X, label="Out: WG_X (Signal)")
axs[0, 0].plot(wavelengths_nm, F_TE_to_WG_Y, label="Out: WG_Y (Crosstalk)")
axs[0, 0].set_title('1. Input: Fiber (TE / X-Pol)')
axs[0, 0].set_ylabel('Gain (dB)')
axs[0, 0].legend()
axs[0, 0].grid(True)

# Plot 2: Forward TM
axs[0, 1].plot(wavelengths_nm, F_TM_to_WG_Y, label="Out: WG_Y (Signal)")
axs[0, 1].plot(wavelengths_nm, F_TM_to_WG_X, label="Out: WG_X (Crosstalk)")
axs[0, 1].set_title('2. Input: Fiber (TM / Y-Pol)')
axs[0, 1].legend()
axs[0, 1].grid(True)

# Plot 3: Reverse WG_X
axs[1, 0].plot(wavelengths_nm, WG_X_to_F_TE, label="Out: Fiber TE (Signal)")
axs[1, 0].plot(wavelengths_nm, WG_X_to_F_TM, label="Out: Fiber TM (Crosstalk)")
axs[1, 0].set_title('3. Input: WG_X')
axs[1, 0].set_xlabel('Wavelength (nm)')
axs[1, 0].set_ylabel('Gain (dB)')
axs[1, 0].legend()
axs[1, 0].grid(True)

# Plot 4: Reverse WG_Y
axs[1, 1].plot(wavelengths_nm, WG_Y_to_F_TM, label="Out: Fiber TM (Signal)")
axs[1, 1].plot(wavelengths_nm, WG_Y_to_F_TE, label="Out: Fiber TE (Crosstalk)")
axs[1, 1].set_title('4. Input: WG_Y')
axs[1, 1].set_xlabel('Wavelength (nm)')
axs[1, 1].legend()
axs[1, 1].grid(True)

plt.tight_layout()
plt.show()