package gui;

import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.text.AbstractDocument;
import java.awt.*;
import java.awt.event.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.List;
import java.util.ArrayList;

/**
 * MainGUIPanel is the primary GUI panel for configuring and running a simulation
 * using various parameters and options. It includes interactive fields,
 * checkboxes, and control buttons for starting and stopping simulations.
 *
 * This panel integrates with a Python backend and includes basic validation
 * and real-time updates.
 */
public class MainGUIPanel extends JPanel {

    MusicPlayer musicPlayer = new MusicPlayer();

    private JButton _startButton;
    private JButton _stopButton;

    private JTextField[] _textFields;
    private double[] _parameters;
    private int _paramNum = 5;
    private final double[] _minValues = { 100, 1, 0.01, 1, 1 };
    private final double[] _maxValues = { 2000, 1000, 1e10, 10000, 20 };


    private JCheckBox[] _checkBoxes;
    private boolean[] _flags;
    private int _flagNum = 3;

    private PythonListener pythonListener;
    private Thread _simulationThread;

    /**
     * Constructs the main GUI panel and initializes all GUI components.
     * Sets up the layout, input fields, checkboxes, and button actions.
     */
    public MainGUIPanel() {
        setLayout(new GridLayout(10, 2, 5, 5)); // 5 sor (4 mező + gomb), 2 oszlop

        //String[] paramLabels = {"Wavelength (nm):", "Waist (um)", "Intensity (W/cm^2):", "Impulse length (fs):", "Spatial grid length:", "Time grid length:"};
        String[] paramNames = {"Wavelength", "Waist", "Intensity", "Impulse length", "Time grid length"};
        String[] paramUnits = {"(nm)", "(um)", "(W/cm²)", "(fs)", ""};
        String[] checkLabels = {"Convergence test:", "Create plots:", "C.A.T (recommended):"};

        _textFields = new JTextField[_paramNum];
        _parameters = new double[_paramNum];

        _checkBoxes = new JCheckBox[_flagNum];
        _flags = new boolean[_flagNum];


        for (int i = 0; i < _paramNum; i++) {
            String labelText = "<html>" + paramNames[i] + " <font color='gray'>" + paramUnits[i] + "</font></html>";
            add(new JLabel(labelText));

            _textFields[i] = new JTextField();
            ((AbstractDocument) _textFields[i].getDocument()).setDocumentFilter(new NumericFilter(false));
            int index = i;

            // Add listener for automatic update
            _textFields[i].getDocument().addDocumentListener(new DocumentListener() {
                public void insertUpdate(DocumentEvent e) { update(); }
                public void removeUpdate(DocumentEvent e) { update(); }
                public void changedUpdate(DocumentEvent e) { update(); }

                /**
                 * Updates the parameter array when text input changes.
                 * Invalid input is silently ignored.
                 */
                private void update() {
                    try {
                        String inputText = _textFields[index].getText().replace(",", ".");
                        _parameters[index] = Double.parseDouble(inputText);
                    } catch (NumberFormatException ex) {
                        // Invalid input, do not update parameter
                    }
                }
            });
            add(_textFields[i]);
        }

        // Create checkboxes and listeners
        for (int i = 0; i < _flagNum; i++) {
            add(new JLabel(checkLabels[i]));
            _checkBoxes[i] = new JCheckBox();
            int index = i;

            // Update flag value when checkbox state changes
            _checkBoxes[i].addItemListener(e -> {
                _flags[index] = _checkBoxes[index].isSelected();
            });

            add(_checkBoxes[i]);
        }

        // Start and Stop buttons
        _startButton = new JButton("Start");
        add(_startButton);

        _stopButton = new JButton("Stop");
        add(_stopButton);

        // Attach action listeners to buttons
        _startButton.addActionListener(new StartButtonListener());
        _stopButton.addActionListener(new StopButtonListener());
    }

    /**
     * Returns the simulation parameters entered by the user.
     *
     * @return a copy of the parameters array
     */
    public double[] getParameters() {
        return _parameters.clone();
    }

    /**
     * Returns the selected feature flags as boolean array.
     *
     * @return a copy of the flags array
     */
    public boolean[] getFlags() {
        return _flags.clone();
    }

    /**
     * Displays an error dialog with the given message and title.
     *
     * @param message the error message
     * @param title   the title of the error dialog
     */
    private void showError(String message, String title) {
        JOptionPane.showMessageDialog(MainGUIPanel.this, message, title, JOptionPane.ERROR_MESSAGE);
    }


    /**
     * Shows an image in a modal dialog with an option to save it.
     *
     * @param imagePath   the path to the image file
     * @param windowTitle the title of the window
     */
    private void showImageWithSaveOption(String imagePath, String windowTitle) {
        File imageFile = new File(imagePath);
        if (!imageFile.exists()) {
            System.out.println("Image not found: " + imagePath);
            return;
        }

        ImageIcon icon = new ImageIcon(imagePath);
        icon.getImage().flush();
        JLabel imageLabel = new JLabel(icon);

        // Create Close and Save buttons
        JButton okButton = new JButton("Close");
        JButton saveButton = new JButton("Save");

        JPanel buttonPanel = new JPanel();
        buttonPanel.add(okButton);
        buttonPanel.add(saveButton);

        // Create main panel with image and buttons
        JPanel mainPanel = new JPanel(new BorderLayout());
        mainPanel.add(imageLabel, BorderLayout.CENTER);
        mainPanel.add(buttonPanel, BorderLayout.SOUTH);

        // Set up and display the dialog
        JDialog dialog = new JDialog((Frame) null, windowTitle, true);
        dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        dialog.setContentPane(mainPanel);
        dialog.pack();
        dialog.setLocationRelativeTo(null);

        // OK button closes the dialog
        okButton.addActionListener(e -> dialog.dispose());

        // Save button opens file save dialog
        saveButton.addActionListener(e -> {
            JFileChooser fileChooser = new JFileChooser();
            fileChooser.setDialogTitle("Save plot as...");
            fileChooser.setSelectedFile(new File("plot.png")); // default file name
            int userSelection = fileChooser.showSaveDialog(dialog);

            if (userSelection == JFileChooser.APPROVE_OPTION) {
                File destFile = fileChooser.getSelectedFile();
                try {
                    Files.copy(imageFile.toPath(), destFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
                    JOptionPane.showMessageDialog(dialog, "Plot saved to:\n" + destFile.getAbsolutePath());
                } catch (IOException ex) {
                    JOptionPane.showMessageDialog(dialog, "Failed to save file:\n" + ex.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
                }
            }
        });

        dialog.setVisible(true);
    }

    /**
     * Displays an image dialog only if the file exists.
     *
     * @param filepath path to the image file
     * @param title    title of the dialog
     */
    private void showIfExists(String filepath, String title) {
        File file = new File(filepath);
        if (file.exists()) {
            showImageWithSaveOption(filepath, title);
        } else {
            System.out.println("File not found, skipping: " + filepath);
        }
    }

    /**
     * Validates if all input fields are filled and contain valid numbers.
     *
     * @return true if valid, false otherwise
     */
    private boolean validateInputs() {
        for (int i = 0; i < _textFields.length; i++) {
            String text = _textFields[i].getText().trim();
            if (text.isEmpty()) {
                showError("Please fill in all fields!", "Missing data");
                return false;
            }
            try {
                _parameters[i] = Double.parseDouble(text.replace(",", "."));
            } catch (NumberFormatException ex) {
                showError("Invalid number in field " + (i + 1) + "!", "Format error");
                return false;
            }
        }
        return true;
    }

    /**
     * Checks if the input parameters are within predefined acceptable ranges.
     *
     * @return true if all values are within range, false otherwise
     */
    private boolean checkParameterRanges() {
        for (int i = 0; i < _parameters.length; i++) {
            if (_parameters[i] < _minValues[i] || _parameters[i] > _maxValues[i]) {
                showError("Value out of range in field " + (i + 1) +
                        "!\nAllowed range: " + _minValues[i] + " – " + _maxValues[i], "Range error");
                return false;
            }
        }
        return true;
    }

    /**
     * Creates a non-blocking dialog to indicate the simulation is in progress.
     * If the CAT flag is enabled, shows an animated GIF and message.
     *
     * @param catFlag whether to show the CAT animation
     * @return the constructed dialog
     */
    private JDialog createProcessingDialog(boolean catFlag) {
        JDialog dialog = new JDialog((JFrame) SwingUtilities.getWindowAncestor(MainGUIPanel.this), "Processing...", false);
        dialog.setLayout(new BorderLayout());

        if (catFlag) {
            ImageIcon gifIcon = new ImageIcon("dancing_cat.gif");
            JLabel gifLabel = new JLabel(gifIcon);
            JLabel textLabel = new JLabel("The code is running, let's listen to some music...", JLabel.CENTER);

            JPanel panel = new JPanel(new BorderLayout());
            panel.add(gifLabel, BorderLayout.CENTER);
            panel.add(textLabel, BorderLayout.SOUTH);
            dialog.add(panel);
            dialog.setSize(350, 350);
        } else {
            dialog.add(new JLabel("The code is running...", JLabel.CENTER), BorderLayout.CENTER);
            dialog.setSize(250, 100);
        }

        dialog.setLocationRelativeTo(MainGUIPanel.this);
        return dialog;
    }

    /**
     * Runs the simulation logic in a background thread.
     * Depending on the flags, music may play and results may be displayed.
     *
     * @param infoDialog the processing dialog to be dismissed upon completion
     */
    private void runSimulation(JDialog infoDialog) {
        boolean catFlag = _flags[2];
        boolean plotFlag = _flags[1];
        boolean convergenceFlag = _flags[0];

        _simulationThread = new Thread(() -> {
            try {
                if (catFlag) musicPlayer.playMusic("music.wav");
                if (pythonListener != null) {
                    pythonListener.onStartButtonPressed(_parameters, _flags);
                } else {
                    System.err.println("PythonListener not set!");
                }
            } catch (Exception ex) {
                ex.printStackTrace();
            } finally {
                if (catFlag) musicPlayer.stopMusic();

                SwingUtilities.invokeLater(() -> {
                    infoDialog.dispose();
                    if (plotFlag) {
                        showIfExists("results/d_plot.png", "Dipole Time Series");
                        showIfExists("results/fft_plot.png", "Dipole Spectrum");
                    }

                    if(convergenceFlag){
                        showIfExists("results/convergence_plot.png", "Convergence Test");
                    }
                });
            }
        });

        _simulationThread.start();
        infoDialog.setVisible(true);
    }

    private void stopSimulation() {
        if (_simulationThread != null && _simulationThread.isAlive()) {
            _simulationThread.interrupt();

            if (pythonListener != null) {
                pythonListener.requestStop();  // Python-nak is szólunk
            }

            musicPlayer.stopMusic();

            JOptionPane.showMessageDialog(this, "Simulation interrupted by user.");
        }
    }




    /// LISTENER ///

    /**
     * ActionListener for the Start button. Validates input and launches simulation.
     */
    public class StartButtonListener implements ActionListener {
        public void actionPerformed(ActionEvent e) {
            if (!validateInputs()) return;
            if (!checkParameterRanges()) return;

            // Confirmation dialog
            int choice = JOptionPane.showOptionDialog(
                    MainGUIPanel.this,
                    "Please make sure to save any important data from the 'results/' folder.\n" +
                            "Running the simulation will overwrite and delete its contents.",
                    "Warning: Results Overwrite",
                    JOptionPane.YES_NO_OPTION,
                    JOptionPane.WARNING_MESSAGE,
                    null,
                    new String[]{"Understood", "Cancel"},
                    "Cancel"
            );

            // Only proceed if user clicks "Understood"
            if (choice != JOptionPane.YES_OPTION) {
                System.out.println("Simulation aborted by user.");
                return;
            }

            boolean catFlag = _flags[2];
            JDialog dialog = createProcessingDialog(catFlag);
            runSimulation(dialog);
        }
    }

    /**
     * ActionListener for the Stop button. Interrupts simulation thread and stops music.
     */
    public class StopButtonListener implements ActionListener {
        public void actionPerformed(ActionEvent e) {
            stopSimulation();
        }
    }


    public void setParametersFromStrings(String[] values) {
        for (int i = 0; i < Math.min(values.length, _textFields.length); i++) {
            _textFields[i].setText(values[i]);
        }
    }

    /**
     * Sets the listener to be called when the Start button is pressed.
     *
     * @param listener the PythonListener to attach
     */
    public void setPythonListener(PythonListener listener) {
        this.pythonListener = listener;
    }
}
