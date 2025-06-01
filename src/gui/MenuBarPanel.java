package gui;

import javax.swing.*;

/**
 * MenuBarPanel defines the menu bar for the main application window.
 * It provides access to tools for saving and loading simulation parameters.
 */
public class MenuBarPanel extends JMenuBar {

    /**
     * Constructs the menu bar with "Tools" menu containing options
     * to save and load simulation parameters.
     *
     * @param parentFrame    the main application frame used for dialog parenting
     * @param mainGUIPanel   the main GUI panel providing parameter data
     */
    public MenuBarPanel(JFrame parentFrame, MainGUIPanel mainGUIPanel) {
        JMenu toolsMenu = new JMenu("Tools");

        // Menu item for saving parameters
        JMenuItem saveParametersItem = new JMenuItem("Save Parameters");
        saveParametersItem.addActionListener(e -> {
            SaveParametersPanel.openSaveParametersDialog(parentFrame, mainGUIPanel);
        });

        // Menu item for loading parameters
        JMenuItem loadParametersItem = new JMenuItem("Load Parameters");
        loadParametersItem.addActionListener(e -> {
            LoadParametersPanel.openLoadParametersDialog(parentFrame, mainGUIPanel);
        });

        toolsMenu.add(saveParametersItem);
        toolsMenu.add(loadParametersItem);
        add(toolsMenu);
    }
}
