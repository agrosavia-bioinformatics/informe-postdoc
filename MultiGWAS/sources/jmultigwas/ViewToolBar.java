/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jmultigwas;

import java.awt.Component;
import javax.swing.JOptionPane;

/**
 *
 * @author lg
 */
public class ViewToolBar extends javax.swing.JPanel {
    Controller controller;

    public ViewToolBar(Controller controller) {
        this.controller = controller;
        initComponents();
    }
    
    public String getToolsToRun () {
        String tools = "";
        if (cboxGwaspoly.isSelected())
                tools += "GWASpoly ";
        if (cboxShesis.isSelected())
                tools += "SHEsis ";
        if (cboxPlink.isSelected())
                tools += "PLINK ";
        if (cboxTassel.isSelected())
                tools += "TASSEL ";
        return (tools);
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jSeparator7 = new javax.swing.JSeparator();
        jToolBar1 = new javax.swing.JToolBar();
        panelToolBar = new javax.swing.JPanel();
        panelTools = new javax.swing.JPanel();
        panelTetra = new javax.swing.JPanel();
        jLabel1 = new javax.swing.JLabel();
        jPanel5 = new javax.swing.JPanel();
        cboxGwaspoly = new javax.swing.JCheckBox();
        cboxShesis = new javax.swing.JCheckBox();
        panelDiplo = new javax.swing.JPanel();
        jLabel8 = new javax.swing.JLabel();
        jPanel17 = new javax.swing.JPanel();
        cboxPlink = new javax.swing.JCheckBox();
        cboxTassel = new javax.swing.JCheckBox();
        buttonRun = new javax.swing.JButton();
        jButton2 = new javax.swing.JButton();
        jButton3 = new javax.swing.JButton();
        labelMultiGWAS = new javax.swing.JLabel();

        setLayout(new java.awt.CardLayout());

        jToolBar1.setBackground(new java.awt.Color(153, 255, 153));
        jToolBar1.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
        jToolBar1.setFloatable(false);
        jToolBar1.setOrientation(javax.swing.SwingConstants.VERTICAL);
        jToolBar1.setPreferredSize(new java.awt.Dimension(155, 50));

        panelToolBar.setPreferredSize(new java.awt.Dimension(150, 140));
        panelToolBar.setRequestFocusEnabled(false);
        panelToolBar.setLayout(new java.awt.BorderLayout());

        panelTools.setPreferredSize(new java.awt.Dimension(150, 120));
        panelTools.setRequestFocusEnabled(false);

        panelTetra.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
        panelTetra.setLayout(new java.awt.BorderLayout());

        jLabel1.setBackground(java.awt.Color.green);
        jLabel1.setText("Tetraploid tools:");
        jLabel1.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
        jLabel1.setOpaque(true);
        panelTetra.add(jLabel1, java.awt.BorderLayout.PAGE_START);

        jPanel5.setLayout(new java.awt.GridLayout(2, 0));

        cboxGwaspoly.setSelected(true);
        cboxGwaspoly.setText("GWASpoly");
        jPanel5.add(cboxGwaspoly);

        cboxShesis.setSelected(true);
        cboxShesis.setText("SHEsis");
        jPanel5.add(cboxShesis);

        panelTetra.add(jPanel5, java.awt.BorderLayout.CENTER);

        panelDiplo.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
        panelDiplo.setLayout(new java.awt.BorderLayout());

        jLabel8.setBackground(java.awt.Color.green);
        jLabel8.setText("Diploid tools:");
        jLabel8.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
        jLabel8.setOpaque(true);
        panelDiplo.add(jLabel8, java.awt.BorderLayout.PAGE_START);

        jPanel17.setLayout(new java.awt.GridLayout(2, 0));

        cboxPlink.setSelected(true);
        cboxPlink.setText("PLINK");
        jPanel17.add(cboxPlink);

        cboxTassel.setSelected(true);
        cboxTassel.setText("TASSEL");
        jPanel17.add(cboxTassel);

        panelDiplo.add(jPanel17, java.awt.BorderLayout.CENTER);

        buttonRun.setText("Run");
        buttonRun.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                buttonRunActionPerformed(evt);
            }
        });

        jButton2.setText("Open...");
        jButton2.setEnabled(false);

        jButton3.setText("Save...");
        jButton3.setEnabled(false);

        javax.swing.GroupLayout panelToolsLayout = new javax.swing.GroupLayout(panelTools);
        panelTools.setLayout(panelToolsLayout);
        panelToolsLayout.setHorizontalGroup(
            panelToolsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelToolsLayout.createSequentialGroup()
                .addGroup(panelToolsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(panelToolsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                        .addComponent(panelDiplo, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(panelTetra, javax.swing.GroupLayout.DEFAULT_SIZE, 160, Short.MAX_VALUE))
                    .addGroup(panelToolsLayout.createSequentialGroup()
                        .addGap(30, 30, 30)
                        .addGroup(panelToolsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(jButton2, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(buttonRun, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jButton3, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))))
                .addContainerGap(218, Short.MAX_VALUE))
        );
        panelToolsLayout.setVerticalGroup(
            panelToolsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelToolsLayout.createSequentialGroup()
                .addGap(5, 5, 5)
                .addComponent(panelTetra, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(10, 10, 10)
                .addComponent(panelDiplo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(buttonRun, javax.swing.GroupLayout.PREFERRED_SIZE, 25, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jButton2)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jButton3)
                .addContainerGap())
        );

        panelToolBar.add(panelTools, java.awt.BorderLayout.CENTER);

        labelMultiGWAS.setBackground(new java.awt.Color(102, 0, 153));
        labelMultiGWAS.setFont(new java.awt.Font("Dialog", 0, 18)); // NOI18N
        labelMultiGWAS.setForeground(java.awt.Color.yellow);
        labelMultiGWAS.setText("MultiGWAS 1.0");
        labelMultiGWAS.setOpaque(true);
        labelMultiGWAS.setPreferredSize(new java.awt.Dimension(82, 55));
        labelMultiGWAS.setRequestFocusEnabled(false);
        panelToolBar.add(labelMultiGWAS, java.awt.BorderLayout.PAGE_START);

        jToolBar1.add(panelToolBar);

        add(jToolBar1, "card2");
    }// </editor-fold>//GEN-END:initComponents

    private void buttonRunActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_buttonRunActionPerformed
        // TODO add your handling code here:
        controller.onRunApplication();
    }//GEN-LAST:event_buttonRunActionPerformed


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton buttonRun;
    private javax.swing.JCheckBox cboxGwaspoly;
    private javax.swing.JCheckBox cboxPlink;
    private javax.swing.JCheckBox cboxShesis;
    private javax.swing.JCheckBox cboxTassel;
    private javax.swing.JButton jButton2;
    private javax.swing.JButton jButton3;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel8;
    private javax.swing.JPanel jPanel17;
    private javax.swing.JPanel jPanel5;
    private javax.swing.JSeparator jSeparator7;
    private javax.swing.JToolBar jToolBar1;
    private javax.swing.JLabel labelMultiGWAS;
    private javax.swing.JPanel panelDiplo;
    private javax.swing.JPanel panelTetra;
    private javax.swing.JPanel panelToolBar;
    private javax.swing.JPanel panelTools;
    // End of variables declaration//GEN-END:variables
}
