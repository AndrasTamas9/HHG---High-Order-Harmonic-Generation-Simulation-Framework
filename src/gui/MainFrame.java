package gui;

import javax.swing.*;

/**
 * MainFrame is the main application window for the HHG (High Harmonic Generation) application.
 * It initializes and holds the central GUI panel and menu bar.
 */
public class MainFrame extends JFrame {

    // The primary panel handling simulation controls and user input
    private MainGUIPanel mainGUIPanel;  // <-- OSZTÁLYSZINTŰ mező lett

    /**
     * Constructs the main frame, sets up window properties,
     * attaches the main panel and menu bar.
     */
    public MainFrame() {
        setTitle("HHG Application");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setSize(400, 400);

        // Initialize the main GUI panel
        mainGUIPanel = new MainGUIPanel();

        // Set up the menu bar from a separate class
        setJMenuBar(new MenuBarPanel(this, mainGUIPanel));

        // Add the main GUI panel to the content pane
        add(mainGUIPanel);

        setVisible(true);
    }

    /**
     * Provides access to the MainGUIPanel instance for external interaction.
     *
     * @return the MainGUIPanel instance
     */
    public MainGUIPanel getMainGUIPanel() {
        return mainGUIPanel;
    }
}
