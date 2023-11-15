/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.jlab.clas.analysis;

import java.awt.Dimension;
import java.awt.Toolkit;
import javax.swing.JFrame;
import org.jlab.clas.swimtools.MagFieldsEngine;
import org.jlab.logging.DefaultLogger;

/**
 *
 * @author veronique
 */
public class AnalysisEngine extends Analysis {
    
    public AnalysisEngine() {
        super("ANAL");
    }
    
    public static void main(String[] args) {

        DefaultLogger.debug();
        
        MagFieldsEngine enf = new MagFieldsEngine();
        enf.init();
        
        AnalysisEngine ae = new AnalysisEngine();
        ae.loadConfiguration();
        
        
    }

}
