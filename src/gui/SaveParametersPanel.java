package gui;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * SaveParametersPanel provides a GUI interface to save simulation
 * parameters from the MainGUIPanel to an external file.
 *
 * Users can select a folder and export the current parameters in XML format.
 */
public class SaveParametersPanel extends JPanel {

    private static JDialog settingsDialog = null;
    private MainGUIPanel _mainGUIPanel;

    // Names and units of parameters
    String[] paramNames = {"Wavelength", "Waist", "Intensity", "Impulse length", "Time grid length"};
    String[] paramUnits = {"(nm)", "(um)", "(W/cm²)", "(fs)", ""};
    String[] checkLabels = {"Convergence test:", "Create plots:", "CAT:"};

    /**
     * Constructs a SaveParametersPanel linked to the given main GUI panel.
     *
     * @param mainGUIPanel the main GUI panel providing the parameters to save
     */
    public SaveParametersPanel(MainGUIPanel mainGUIPanel) {
        _mainGUIPanel = mainGUIPanel;

        setLayout(new GridLayout(2, 1, 5, 5));

        add(new JLabel("Select a folder to save:"));

        JButton saveButton = new JButton("Save parameters");
        add(saveButton);

        // Button triggers parameter saving
        saveButton.addActionListener(e -> saveParameters());
    }

    /**
     * Opens a file chooser to select save destination and writes parameters to XML.
     */
    private void saveParameters() {
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);

        int result = fileChooser.showSaveDialog(this);

        if (result == JFileChooser.APPROVE_OPTION) {
            File selectedDirectory = fileChooser.getSelectedFile();

            // Default file name with timestamp
            String timestamp = new java.text.SimpleDateFormat("yyyy-MM-dd_HH-mm-ss").format(new java.util.Date());
            String defaultName = "HHG_parameters_" + timestamp;
            String input = JOptionPane.showInputDialog(this, "Enter the file name (without extension):", defaultName);

            String fileName;
            if (input == null || input.trim().isEmpty()) {
                fileName = defaultName + ".xml";
            } else {
                fileName = input.trim() + ".xml";
            }

            String path = selectedDirectory.getAbsolutePath() + File.separator + fileName;

            System.out.println("Saving to:" + path);

            double[] parameters = _mainGUIPanel.getParameters();
            boolean[] flags = _mainGUIPanel.getFlags();

            try (FileWriter writer = new FileWriter(path)) {
                // Begin XML document
                writer.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
                writer.write("<SimulationParameters>\n");
                // Write parameters block
                writer.write("  <Parameters>\n");
                for (int i = 0; i < parameters.length; i++) {
                    String name = paramNames[i].replace(" ", "");
                    String unit = paramUnits[i].replaceAll("[()]", "");  // "nm", "um", stb.
                    writer.write("    <" + name + (unit.isEmpty() ? "" : " unit=\"" + unit + "\"") + ">");
                    writer.write(String.valueOf(parameters[i]));
                    writer.write("</" + name + ">\n");
                }
                writer.write("  </Parameters>\n");

                // Write flags block
                writer.write("  <Flags>\n");
                for (int i = 0; i < flags.length; i++) {
                    String label = checkLabels[i].replace(":", "").replaceAll("\\s+", "");
                    writer.write("    <" + label + ">" + flags[i] + "</" + label + ">\n");
                }
                writer.write("  </Flags>\n");
                writer.write("</SimulationParameters>\n");

                // Notify user of success
                JOptionPane.showMessageDialog(this, "Save successful!\n" + path);
                settingsDialog.dispose();

            } catch (IOException ex) {
                ex.printStackTrace();
                JOptionPane.showMessageDialog(this, "Error saving the file!", "Error", JOptionPane.ERROR_MESSAGE);
            }
        }
    }


    /**
     * Opens the Save Parameters dialog if it is not already open.
     * Ensures only one instance of the dialog is shown at a time.
     *
     * @param parent        the parent window for modal alignment
     * @param mainGUIPanel  the main panel providing parameters to save
     */
    public static void openSaveParametersDialog(JFrame parent, MainGUIPanel mainGUIPanel) {
        if (settingsDialog == null || !settingsDialog.isShowing()) {
            settingsDialog = new JDialog(parent, "Save Parameters", true);
            settingsDialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
            settingsDialog.setSize(400, 150);
            settingsDialog.setLocationRelativeTo(parent);
            settingsDialog.add(new SaveParametersPanel(mainGUIPanel));
            settingsDialog.setVisible(true);
        } else {
            settingsDialog.toFront();
        }
    }
}
