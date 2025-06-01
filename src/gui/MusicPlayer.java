package gui;

import javax.sound.sampled.*;
import java.io.File;
import java.io.IOException;

/**
 * MusicPlayer provides basic functionality to play and stop audio files.
 * It uses the javax.sound.sampled API to stream WAV audio files and
 * can loop playback continuously.
 */
public class MusicPlayer {
    private Clip clip;

    /**
     * Plays the specified audio file in a continuous loop.
     *
     * @param filePath the path to the audio file (must be a supported format like WAV)
     */
    public void playMusic(String filePath) {
        try {
            AudioInputStream audioStream = AudioSystem.getAudioInputStream(new File(filePath));
            clip = AudioSystem.getClip();
            clip.open(audioStream);
            clip.loop(Clip.LOOP_CONTINUOUSLY); // loop playback indefinitely
        } catch (UnsupportedAudioFileException | IOException | LineUnavailableException e) {
            e.printStackTrace();
        }
    }

    /**
     * Stops audio playback and releases system resources.
     */
    public void stopMusic() {
        if (clip != null) {
            clip.stop();
            clip.close();
        }
    }
}
