package gui;

/**
 * PythonListener is an interface for communication between the GUI
 * and a Python backend or simulation controller.
 *
 * Implementations of this interface handle the start and stop
 * requests initiated from the GUI.
 */
public interface PythonListener {
    /**
     * Called when the user presses the Start button in the GUI.
     * Should initiate a computation or simulation using the provided parameters.
     *
     * @param parameters array of simulation parameters (e.g. wavelength, intensity)
     * @param flags boolean array representing feature toggles (e.g. plotting enabled)
     */
    void onStartButtonPressed(double[] parameters, boolean[] flags);

    /**
     * Called when the user requests to stop the ongoing simulation.
     * Should trigger graceful shutdown of the Python-side process if running.
     */
    void requestStop();
}
