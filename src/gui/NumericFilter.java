package gui;

import javax.swing.text.*;

/**
 * NumericFilter restricts input in text components to numeric values only.
 * Optionally allows or disallows negative numbers.
 *
 * Can be used with JTextField via AbstractDocument#setDocumentFilter().
 */
public class NumericFilter extends DocumentFilter {

    private boolean allowNegative;

    /**
     * Constructs a NumericFilter allowing negative numbers by default.
     */
    public NumericFilter() {
        this(true);
    }

    /**
     * Constructs a NumericFilter with an option to allow or disallow negatives.
     *
     * @param allowNegative whether negative numbers should be permitted
     */
    public NumericFilter(boolean allowNegative) {
        this.allowNegative = allowNegative;
    }

    /**
     * Validates inserted text before allowing it into the document.
     */
    @Override
    public void insertString(FilterBypass fb, int offset, String string, AttributeSet attr)
            throws BadLocationException {
        String currentText = fb.getDocument().getText(0, fb.getDocument().getLength());
        String newText = currentText.substring(0, offset) + string + currentText.substring(offset);
        if (isNumericInput(newText)) {
            super.insertString(fb, offset, string, attr);
        }
    }

    /**
     * Validates replaced text before allowing it into the document.
     */
    @Override
    public void replace(FilterBypass fb, int offset, int length, String text, AttributeSet attrs)
            throws BadLocationException {
        String currentText = fb.getDocument().getText(0, fb.getDocument().getLength());
        String newText = currentText.substring(0, offset) + text + currentText.substring(offset + length);
        if (isNumericInput(newText)) {
            super.replace(fb, offset, length, text, attrs);
        }
    }

    /**
     * Validates whether the given string is a valid numeric input
     * according to the filter settings (allows optional negatives and decimals).
     *
     * @param text the string to validate
     * @return true if the input is numeric or empty, false otherwise
     */
    private boolean isNumericInput(String text) {
        if (text.isEmpty()) return true; // engedi törléskor is
        String regex = allowNegative ? "-?\\d*(\\.\\d*)?" : "\\d*(\\.\\d*)?";
        return text.matches(regex);
    }
}
