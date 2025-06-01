package gui;


import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.util.ArrayList;


/**
 * LoadParametersPanel provides a simple GUI interface for loading
 * parameter files into the MainGUIPanel.
 *
 * It includes a button to select and load a file, and parses its
 * content to update simulation parameters.
 */
public class LoadParametersPanel extends JPanel {

    private static JDialog settingsDialog = null;
    private MainGUIPanel _mainGUIPanel;

    /**
     * Constructs the LoadParametersPanel and links it with the main GUI.
     *
     * @param mainGUIPanel the main GUI panel to send the loaded parameters to
     */
    public LoadParametersPanel(MainGUIPanel mainGUIPanel) {
        _mainGUIPanel = mainGUIPanel;

        // Set layout with 2 rows and 1 column
        setLayout(new GridLayout(2, 1, 5, 5));

        add(new JLabel("Select a file:"));

        JButton loadButton = new JButton("Load parameters");
        add(loadButton);

        // Button action to trigger file loading
        loadButton.addActionListener(e -> loadParameters());
    }

    /**
     * Opens a file chooser dialog, reads parameters from the selected file,
     * and updates the main GUI panel accordingly.
     */
    private void loadParameters() {
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);

        int result = fileChooser.showOpenDialog(this);

        if (result == JFileChooser.APPROVE_OPTION) {
            File selectedFile = fileChooser.getSelectedFile();

            try {
                // Prepare to read the XML or custom parameter file
                javax.xml.parsers.DocumentBuilderFactory factory = javax.xml.parsers.DocumentBuilderFactory.newInstance();
                javax.xml.parsers.DocumentBuilder builder = factory.newDocumentBuilder();
                org.w3c.dom.Document doc = builder.parse(selectedFile);

                // Select <Parameters> tag from XML
                org.w3c.dom.NodeList paramNodes = doc.getElementsByTagName("Parameters");
                if (paramNodes.getLength() == 0) {
                    throw new Exception("No <Parameters> section found in XML.");
                }

                org.w3c.dom.Element parametersElement = (org.w3c.dom.Element) paramNodes.item(0);
                org.w3c.dom.NodeList children = parametersElement.getChildNodes();

                // Collect values from all element-type child nodes
                ArrayList<String> values = new ArrayList<>();
                for (int i = 0; i < children.getLength(); i++) {
                    org.w3c.dom.Node node = children.item(i);
                    if (node.getNodeType() == org.w3c.dom.Node.ELEMENT_NODE) {
                        String value = node.getTextContent().trim();
                        values.add(value);
                    }
                }

                // Convert collected values into a String array
                String[] paramArray = values.toArray(new String[0]);

                // Pass parameters to the main GUI panel
                _mainGUIPanel.setParametersFromStrings(paramArray);

                // Notify user and close the dialog
                JOptionPane.showMessageDialog(this, "Parameters loaded!");
                settingsDialog.dispose();

            } catch (Exception ex) {
                ex.printStackTrace();
                JOptionPane.showMessageDialog(this, "Error loading XML!", "Error", JOptionPane.ERROR_MESSAGE);
            }
        }
    }



    /**
     * Opens the Load Parameters dialog window if not already open.
     * This is a static method to ensure only one instance exists at a time.
     *
     * @param parent         the parent frame for modal dialog alignment
     * @param mainGUIPanel   the main GUI panel to which parameters will be applied
     */
    public static void openLoadParametersDialog(JFrame parent, MainGUIPanel mainGUIPanel) {
        if (settingsDialog == null || !settingsDialog.isShowing()) {
            settingsDialog = new JDialog(parent, "Load Parameters", true);
            settingsDialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
            settingsDialog.setSize(400, 150);
            settingsDialog.setLocationRelativeTo(parent);

            // Embed the panel into the dialog
            settingsDialog.add(new LoadParametersPanel(mainGUIPanel));
            settingsDialog.setVisible(true);
        } else {
            settingsDialog.toFront();
        }
    }
}
