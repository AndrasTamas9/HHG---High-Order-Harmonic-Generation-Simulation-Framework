# ─── IMPORTS ────────────────────────────────────────────────────────────────
# ─── IMPORTS ────────────────────────────────────────────────────────────────
import jpype
import jpype.imports
from jpype import JProxy, java
import os
import time
import matplotlib
matplotlib.use('Agg')
from hhg_python import run_simulation, cleanup_result_images, cleanup_result_files, run_convergence_test
from hhg_python.interface import hhg_lib, convergence_lib


# ─── JAR AND JVM ───────────────────────────────────────────────────────
jar_path = os.path.abspath("gui.jar")

if not jpype.isJVMStarted():
    jpype.startJVM(classpath=[jar_path])

from gui import MainFrame, PythonListener


# ─── PYTHON - JAVA Listener ─────────────────────
class Listener:
    def __init__(self):
        self.stop_requested = False

    def onStartButtonPressed(self, parameters, flags):
        self.stop_requested = False
        
        cleanup_result_files()

        print(">>> Java meghívta a Python oldalt")
        print("Paraméterek:", list(parameters))
        print("Flag-ek:", list(flags))
        print("Összeg:", sum(parameters))
        
        try:
            os.remove("output_success.txt")
        except FileNotFoundError:
            pass
            
        # ⬇️ Paraméterek kicsomagolása
        wavelength = float(parameters[0])   # nm
        waist = float(parameters[1])       # μm
        intensity = float(parameters[2])   # x10^14 W/cm^2
        T = float(parameters[3])           # fs
        exponent = int(parameters[4])      # 2^exponent = TSIZE

        print(">>> Python: szimuláció indítása a C++ kóddal...")
        
        if(flags[0] == False):
            run_simulation(
                wavelength=wavelength,
                waist=waist,
                intensity=intensity,
                pulse_length=T,
                exponent=exponent,
                generate_plots=flags[1]
            )
        else:
            run_convergence_test(
                wavelength=wavelength,
                waist=waist,
                intensity=intensity,
                pulse_length=T,
                exponent=exponent,
                generate_plots=flags[1]
            )
            

        print(">>> Python: C++ szimuláció befejeződött.")

    def requestStop(self):
        print(">>> Python: stop jelet kapott")
        self.stop_requested = True
        hhg_lib.requestStop()
        convergence_lib.requestStop()

        

# ─── START GUI ──────────────────────────────────────────────────────────
proxy = JProxy(PythonListener, inst=Listener())

frame = MainFrame()
panel = frame.getMainGUIPanel()
panel.setPythonListener(proxy)

frame.setVisible(True)



# ─── AFTER CLOSING THE GUI ────────────────────────────────────────────
def on_window_close():
    print(">>> GUI closed – cleaning up...")
    cleanup_result_images()

def on_window_closed():
    print(">>> GUI fully closed – shutting down JVM.")
    jpype.shutdownJVM()

frame.addWindowListener(
    JProxy(
        java.awt.event.WindowListener,
        inst=type("WindowCloser", (), {
            "windowClosing": lambda self, e: on_window_close(),
            "windowClosed": lambda self, e: on_window_closed(),
            "windowOpened": lambda self, e: None,
            "windowIconified": lambda self, e: None,
            "windowDeiconified": lambda self, e: None,
            "windowActivated": lambda self, e: None,
            "windowDeactivated": lambda self, e: None,
        })()
    )
)


