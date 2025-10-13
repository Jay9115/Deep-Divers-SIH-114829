// Analysis Dashboard JavaScript functionality

document.addEventListener('DOMContentLoaded', function() {
    // Initialize analysis dashboard features
    initAnalysisDashboard();
    initZoneGrid();
    initAnalysisProgress();
    initRealTimeUpdates();
    initControls();
});

// Initialize main dashboard functionality
function initAnalysisDashboard() {
    // Animate main statistics
    const statNumbers = document.querySelectorAll('.stat-number-large');
    
    const animateCounter = (element) => {
        const target = parseInt(element.getAttribute('data-target'));
        const duration = 2000;
        const step = target / (duration / 16);
        let current = 0;

        const timer = setInterval(() => {
            current += step;
            if (current >= target) {
                current = target;
                clearInterval(timer);
            }
            
            let displayValue = Math.floor(current);
            if (target > 1000) {
                displayValue = displayValue.toLocaleString();
            }
            
            // Add % for percentage values
            const label = element.nextElementSibling.textContent;
            if (label.includes('Confidence') || label.includes('Avg')) {
                displayValue += '%';
            }
            
            element.textContent = displayValue;
        }, 16);
    };

    // Intersection Observer for stat animation
    const observer = new IntersectionObserver((entries) => {
        entries.forEach(entry => {
            if (entry.isIntersecting && !entry.target.hasAttribute('data-animated')) {
                entry.target.setAttribute('data-animated', 'true');
                animateCounter(entry.target);
            }
        });
    });

    statNumbers.forEach(stat => observer.observe(stat));
}

// Initialize zone grid system
function initZoneGrid() {
    const zoneGrid = document.getElementById('zoneGrid');
    const totalZones = 64;
    
    // Generate zone grid (8x8)
    for (let i = 1; i <= totalZones; i++) {
        const zoneCell = document.createElement('div');
        zoneCell.className = 'zone-cell';
        zoneCell.textContent = i;
        zoneCell.setAttribute('data-zone', i);
        
        // Add click handler
        zoneCell.addEventListener('click', () => toggleZone(i));
        
        zoneGrid.appendChild(zoneCell);
    }
    
    // Simulate active zones (12, 34, 7 as shown in images)
    setTimeout(() => {
        activateZone(12);
        activateZone(34);
        activateZone(7);
    }, 1000);
}

// Toggle zone selection
function toggleZone(zoneNumber) {
    const zoneCell = document.querySelector(`[data-zone="${zoneNumber}"]`);
    
    if (zoneCell.classList.contains('active')) {
        zoneCell.classList.remove('active');
    } else {
        zoneCell.classList.add('active');
        // Simulate search start
        setTimeout(() => {
            zoneCell.classList.add('searching');
        }, 500);
    }
}

// Activate specific zone
function activateZone(zoneNumber) {
    const zoneCell = document.querySelector(`[data-zone="${zoneNumber}"]`);
    if (zoneCell) {
        zoneCell.classList.add('active');
        
        // Start searching animation
        setTimeout(() => {
            zoneCell.classList.add('searching');
        }, 500);
        
        // Complete search after some time
        setTimeout(() => {
            zoneCell.classList.remove('searching');
            zoneCell.classList.add('completed');
        }, 3000 + Math.random() * 2000);
    }
}

// Initialize analysis progress
function initAnalysisProgress() {
    const steps = document.querySelectorAll('.step');
    let currentStep = 3; // Currently on step 4 (0-indexed)
    
    // Simulate progress through steps
    setTimeout(() => {
        progressToNextStep(currentStep + 1);
    }, 3000);
    
    setTimeout(() => {
        progressToNextStep(currentStep + 2);
    }, 6000);
    
    setTimeout(() => {
        progressToNextStep(currentStep + 3);
        showAnalysisComplete();
    }, 9000);
}

// Progress to next step
function progressToNextStep(stepIndex) {
    const steps = document.querySelectorAll('.step');
    if (steps[stepIndex]) {
        steps[stepIndex].classList.add('active');
        
        // Mark previous step as completed
        if (steps[stepIndex - 1]) {
            steps[stepIndex - 1].classList.add('completed');
        }
        
        // Update step badge
        const stepBadge = document.querySelector('.step-badge');
        if (stepBadge) {
            stepBadge.textContent = `Processing step ${stepIndex + 1}...`;
        }
    }
}

// Show analysis complete
function showAnalysisComplete() {
    const stepBadge = document.querySelector('.step-badge');
    if (stepBadge) {
        stepBadge.textContent = 'Analysis Complete!';
        stepBadge.style.background = 'var(--gradient-secondary)';
    }
    
    // Mark all steps as completed
    const steps = document.querySelectorAll('.step');
    steps.forEach(step => {
        step.classList.add('completed');
    });
    
    // Show output layer
    const outputLayer = document.getElementById('outputLayer');
    if (outputLayer) {
        setTimeout(() => {
            outputLayer.style.display = 'block';
            outputLayer.scrollIntoView({ behavior: 'smooth', block: 'start' });
        }, 1000);
    }
    
    // Show completion notification
    showNotification('Analysis completed successfully! Results are ready for download.', 'success');
}

// Real-time updates for sequence processing
function initRealTimeUpdates() {
    const sequences = document.querySelectorAll('.sequence-item');
    
    // Simulate real-time sequence processing
    setInterval(() => {
        updateRandomSequence();
    }, 2000);
    
    // Update zone progress
    setInterval(() => {
        updateZoneProgress();
    }, 1500);
}

// Update random sequence status
function updateRandomSequence() {
    const processingSequences = document.querySelectorAll('.sequence-item.processing');
    
    if (processingSequences.length > 0) {
        const randomSeq = processingSequences[Math.floor(Math.random() * processingSequences.length)];
        const matchStatus = randomSeq.querySelector('.match-status');
        
        // Generate random match percentage
        const matchPercent = Math.floor(Math.random() * 100);
        matchStatus.textContent = `${matchPercent}% match`;
        
        // Update class based on match quality
        randomSeq.classList.remove('processing');
        if (matchPercent >= 70) {
            randomSeq.classList.add('match');
        } else if (matchPercent >= 50) {
            randomSeq.classList.add('match');
        } else {
            randomSeq.classList.add('low-match');
        }
    }
}

// Update zone progress
function updateZoneProgress() {
    const progressNumbers = document.querySelectorAll('.progress-numbers');
    
    progressNumbers.forEach(progress => {
        const text = progress.textContent;
        const match = text.match(/(\d+) \/ (\d+) sequences (\d+)% complete/);
        
        if (match) {
            let current = parseInt(match[1]);
            const total = parseInt(match[2]);
            let percent = parseInt(match[3]);
            
            if (current < total) {
                current = Math.min(current + Math.floor(Math.random() * 20), total);
                percent = Math.floor((current / total) * 100);
                progress.textContent = `${current} / ${total} sequences ${percent}% complete`;
            }
        }
    });
}

// Initialize control buttons
function initControls() {
    const pauseBtn = document.getElementById('pauseBtn');
    const resetBtn = document.getElementById('resetBtn');
    const zoomInBtn = document.getElementById('zoomInBtn');
    const zoomOutBtn = document.getElementById('zoomOutBtn');
    
    let isPaused = false;
    
    // Pause/Resume functionality
    if (pauseBtn) {
        pauseBtn.addEventListener('click', () => {
            isPaused = !isPaused;
            
            if (isPaused) {
                pauseBtn.innerHTML = '<i class="fas fa-play"></i> Resume Analysis';
                pauseBtn.classList.remove('btn-danger');
                pauseBtn.classList.add('btn-primary');
                
                // Pause all animations
                pauseAnimations();
                showNotification('Analysis paused', 'info');
            } else {
                pauseBtn.innerHTML = '<i class="fas fa-pause"></i> Pause Analysis';
                pauseBtn.classList.remove('btn-primary');
                pauseBtn.classList.add('btn-danger');
                
                // Resume animations
                resumeAnimations();
                showNotification('Analysis resumed', 'success');
            }
        });
    }
    
    // Reset functionality
    if (resetBtn) {
        resetBtn.addEventListener('click', () => {
            if (confirm('Are you sure you want to reset the analysis? All progress will be lost.')) {
                resetAnalysis();
                showNotification('Analysis reset', 'info');
            }
        });
    }
    
    // Zoom controls for zone grid
    if (zoomInBtn && zoomOutBtn) {
        let zoomLevel = 1;
        const zoneGrid = document.getElementById('zoneGrid');
        
        zoomInBtn.addEventListener('click', () => {
            if (zoomLevel < 2) {
                zoomLevel += 0.2;
                zoneGrid.style.transform = `scale(${zoomLevel})`;
                zoneGrid.style.transformOrigin = 'center';
            }
        });
        
        zoomOutBtn.addEventListener('click', () => {
            if (zoomLevel > 0.6) {
                zoomLevel -= 0.2;
                zoneGrid.style.transform = `scale(${zoomLevel})`;
                zoneGrid.style.transformOrigin = 'center';
            }
        });
    }
}

// Pause animations
function pauseAnimations() {
    const processingElements = document.querySelectorAll('.processing, .searching');
    processingElements.forEach(el => {
        el.style.animationPlayState = 'paused';
    });
}

// Resume animations
function resumeAnimations() {
    const processingElements = document.querySelectorAll('.processing, .searching');
    processingElements.forEach(el => {
        el.style.animationPlayState = 'running';
    });
}

// Reset analysis
function resetAnalysis() {
    // Reset steps
    const steps = document.querySelectorAll('.step');
    steps.forEach((step, index) => {
        step.classList.remove('active', 'completed');
        if (index === 0) {
            step.classList.add('active');
        }
    });
    
    // Reset zones
    const zoneCells = document.querySelectorAll('.zone-cell');
    zoneCells.forEach(cell => {
        cell.classList.remove('active', 'searching', 'completed');
    });
    
    // Reset sequences
    const sequences = document.querySelectorAll('.sequence-item');
    sequences.forEach(seq => {
        seq.classList.remove('match', 'low-match');
        seq.classList.add('processing');
        seq.querySelector('.match-status').textContent = 'processing...';
    });
    
    // Reset progress
    const progressNumbers = document.querySelectorAll('.progress-numbers');
    progressNumbers.forEach(progress => {
        const text = progress.textContent;
        const match = text.match(/(\d+) \/ (\d+) sequences/);
        if (match) {
            const total = match[2];
            progress.textContent = `0 / ${total} sequences 0% complete`;
        }
    });
    
    // Reset step badge
    const stepBadge = document.querySelector('.step-badge');
    if (stepBadge) {
        stepBadge.textContent = 'Processing step 1...';
        stepBadge.style.background = 'var(--gradient-primary)';
    }
    
    // Restart analysis
    setTimeout(() => {
        initAnalysisProgress();
        activateZone(12);
        activateZone(34);
        activateZone(7);
    }, 1000);
}

// Enhanced notification system with analysis-specific styling
function showNotification(message, type = 'info') {
    const notification = document.createElement('div');
    notification.className = `notification notification-${type}`;
    notification.innerHTML = `
        <i class="fas fa-${getNotificationIcon(type)}"></i>
        <div class="notification-content">
            <div class="notification-title">${getNotificationTitle(type)}</div>
            <div class="notification-message">${message}</div>
        </div>
        <button class="notification-close">&times;</button>
    `;

    document.body.appendChild(notification);
    
    // Enhanced notification styles
    const style = document.createElement('style');
    style.textContent = `
        .notification {
            position: fixed;
            top: 100px;
            right: 20px;
            background: var(--secondary-dark);
            border: 1px solid ${getNotificationColor(type)};
            border-radius: var(--border-radius);
            padding: 1rem 1.5rem;
            display: flex;
            align-items: flex-start;
            gap: 1rem;
            z-index: 3000;
            animation: slideInRight 0.3s ease-out;
            max-width: 400px;
            box-shadow: var(--shadow-card);
            backdrop-filter: blur(10px);
        }
        
        .notification i {
            color: ${getNotificationColor(type)};
            font-size: 1.5rem;
            margin-top: 0.25rem;
        }
        
        .notification-content {
            flex: 1;
        }
        
        .notification-title {
            font-weight: 600;
            color: var(--text-white);
            margin-bottom: 0.25rem;
        }
        
        .notification-message {
            color: var(--text-gray);
            font-size: 0.9rem;
        }
        
        .notification-close {
            background: none;
            border: none;
            color: var(--text-gray);
            font-size: 1.2rem;
            cursor: pointer;
            padding: 0;
        }
    `;
    
    if (!document.head.querySelector('style[data-analysis-notifications]')) {
        style.setAttribute('data-analysis-notifications', 'true');
        document.head.appendChild(style);
    }

    // Close notification
    notification.querySelector('.notification-close').addEventListener('click', () => {
        notification.remove();
    });

    // Auto remove
    setTimeout(() => {
        if (notification.parentNode) {
            notification.remove();
        }
    }, 5000);
}

// Helper functions for notifications
function getNotificationIcon(type) {
    const icons = {
        success: 'check-circle',
        error: 'exclamation-circle',
        info: 'info-circle',
        warning: 'exclamation-triangle'
    };
    return icons[type] || 'info-circle';
}

function getNotificationTitle(type) {
    const titles = {
        success: 'Success',
        error: 'Error',
        info: 'Information',
        warning: 'Warning'
    };
    return titles[type] || 'Notification';
}

function getNotificationColor(type) {
    const colors = {
        success: 'var(--accent-green)',
        error: 'var(--accent-orange)',
        info: 'var(--accent-blue)',
        warning: 'var(--accent-orange)'
    };
    return colors[type] || 'var(--accent-blue)';
}

// Output Layer Functions
function downloadFile(filename) {
    showNotification(`Downloading ${filename}...`, 'info');
    
    // Simulate file download
    const link = document.createElement('a');
    link.href = '#';
    link.download = filename;
    
    // In a real implementation, this would be the actual file URL
    console.log(`Downloading file: ${filename}`);
    
    setTimeout(() => {
        showNotification(`${filename} downloaded successfully!`, 'success');
    }, 1500);
}

function previewFile(fileType) {
    showNotification(`Opening ${fileType} preview...`, 'info');
    
    // Create preview modal content based on file type
    let previewContent = '';
    
    switch(fileType) {
        case 'asv_table':
            previewContent = `
                <h3>ASV Table Preview</h3>
                <div class="file-preview">
                    <pre>
# Constructed from biom file
#OTU ID	Sample1	Sample2	Sample3	taxonomy
ASV_001	1245	2341	1876	k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria
ASV_002	876	1234	2145	k__Bacteria;p__Bacteroidetes;c__Bacteroidia
ASV_003	543	876	1234	k__Archaea;p__Euryarchaeota;c__Methanomicrobia
...
                    </pre>
                </div>
            `;
            break;
        case 'taxonomy':
            previewContent = `
                <h3>Taxonomy Classification Preview</h3>
                <div class="file-preview">
                    <pre>
Feature ID	Taxon	Confidence
ASV_001	k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Oceanospirillales	0.97
ASV_002	k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales	0.94
ASV_003	k__Archaea; p__Euryarchaeota; c__Methanomicrobia; o__Methanomicrobiales	0.89
...
                    </pre>
                </div>
            `;
            break;
        case 'novel_taxa':
            previewContent = `
                <h3>Novel Taxa Report Preview</h3>
                <div class="file-preview">
                    <h4>Priority Novel Taxa Discoveries</h4>
                    <p><strong>ASV_047:</strong> Potential new species in family Oceanospirillaceae</p>
                    <p>Confidence: Low (0.23) - Requires further investigation</p>
                    <p>Depth: 2,340m - Abyssal zone specimen</p>
                    <hr>
                    <p><strong>ASV_123:</strong> Unclassified deep-sea bacterium</p>
                    <p>Confidence: Very Low (0.12) - Novel discovery candidate</p>
                    <p>Depth: 3,890m - Hadal zone specimen</p>
                </div>
            `;
            break;
        default:
            previewContent = `
                <h3>File Preview</h3>
                <p>Preview for ${fileType} would be displayed here in a real implementation.</p>
            `;
    }
    
    // Show preview in a modal (simplified version)
    const existingModal = document.getElementById('previewModal');
    if (existingModal) {
        existingModal.remove();
    }
    
    const previewModal = document.createElement('div');
    previewModal.id = 'previewModal';
    previewModal.className = 'modal';
    previewModal.style.display = 'block';
    previewModal.innerHTML = `
        <div class="modal-content" style="max-width: 800px;">
            <div class="modal-header">
                <h3>File Preview</h3>
                <span class="modal-close" onclick="document.getElementById('previewModal').remove()">&times;</span>
            </div>
            <div class="modal-body">
                ${previewContent}
            </div>
        </div>
    `;
    
    document.body.appendChild(previewModal);
}

function openVisualization(type) {
    showNotification('Opening interactive visualization...', 'info');
    // In a real implementation, this would open the actual visualization
    window.open('#', '_blank');
}

function openREnvironment() {
    showNotification('Preparing R environment...', 'info');
    // In a real implementation, this would integrate with RStudio or similar
    console.log('Opening R environment with phyloseq data');
}

function viewLog() {
    showNotification('Opening processing log...', 'info');
    
    const logContent = `
        <h3>Processing Log</h3>
        <div class="file-preview">
            <pre>
[2025-10-07 14:23:15] Analysis started
[2025-10-07 14:23:16] Input file: uploaded_seq_001.fasta (268 bp)
[2025-10-07 14:23:17] ZHNSW algorithm initialized
[2025-10-07 14:23:18] Zone partitioning: 64 zones created
[2025-10-07 14:23:20] Quality assessment: PASSED
[2025-10-07 14:23:22] Taxonomic classification: 98.7% confidence
[2025-10-07 14:23:24] Novel species detection: 2 candidates found
[2025-10-07 14:23:25] Analysis completed successfully
            </pre>
        </div>
    `;
    
    const logModal = document.createElement('div');
    logModal.className = 'modal';
    logModal.style.display = 'block';
    logModal.innerHTML = `
        <div class="modal-content">
            <div class="modal-header">
                <h3>Processing Log</h3>
                <span class="modal-close" onclick="this.closest('.modal').remove()">&times;</span>
            </div>
            <div class="modal-body">
                ${logContent}
            </div>
        </div>
    `;
    
    document.body.appendChild(logModal);
}

function downloadFigures() {
    showNotification('Preparing figure package...', 'info');
    setTimeout(() => {
        showNotification('Figure package downloaded successfully!', 'success');
    }, 2000);
}

function previewVisualization(type) {
    showNotification('Loading interactive visualization preview...', 'info');
    
    const vizContent = `
        <h3>Interactive Dashboard Preview</h3>
        <div class="viz-preview-container">
            <div class="viz-preview-grid">
                <div class="viz-preview-item">
                    <i class="fas fa-project-diagram"></i>
                    <h4>UMAP Embeddings</h4>
                    <p>2D/3D dimensional reduction visualization</p>
                </div>
                <div class="viz-preview-item">
                    <i class="fas fa-chart-area"></i>
                    <h4>Diversity Plots</h4>
                    <p>Alpha & Beta diversity analysis</p>
                </div>
                <div class="viz-preview-item">
                    <i class="fas fa-sun"></i>
                    <h4>Taxonomy Sunburst</h4>
                    <p>Hierarchical taxonomic representation</p>
                </div>
                <div class="viz-preview-item">
                    <i class="fas fa-chart-bar"></i>
                    <h4>Taxonomy Stats</h4>
                    <p>Statistical summaries and metrics</p>
                </div>
            </div>
            <div style="text-align: center; padding: 1rem; background: rgba(0, 255, 136, 0.1); border-radius: 8px; border: 1px solid rgba(0, 255, 136, 0.3);">
                <p style="color: var(--accent-green); margin: 0;"><i class="fas fa-info-circle"></i> 
                Interactive dashboard with real-time filtering, zooming, and data exploration capabilities</p>
            </div>
        </div>
    `;
    
    const vizModal = document.createElement('div');
    vizModal.className = 'modal';
    vizModal.style.display = 'block';
    vizModal.innerHTML = `
        <div class="modal-content" style="max-width: 800px;">
            <div class="modal-header">
                <h3>Interactive Visualization Preview</h3>
                <span class="modal-close" onclick="this.closest('.modal').remove()">&times;</span>
            </div>
            <div class="modal-body">
                ${vizContent}
            </div>
        </div>
    `;
    
    document.body.appendChild(vizModal);
}

function previewFigures() {
    showNotification('Loading figure gallery...', 'info');
    
    const figureContent = `
        <h3>Publication-Ready Figures Gallery</h3>
        <div class="figure-gallery">
            <div class="figure-item">
                <h4>UMAP Plot (300 DPI)</h4>
                <div class="figure-placeholder">
                    <i class="fas fa-project-diagram" style="font-size: 3rem; color: var(--accent-blue);"></i>
                    <p>High-resolution UMAP embedding visualization</p>
                </div>
            </div>
            <div class="figure-item">
                <h4>Taxonomy Sunburst (PDF)</h4>
                <div class="figure-placeholder">
                    <i class="fas fa-sun" style="font-size: 3rem; color: var(--accent-orange);"></i>
                    <p>Vector-based taxonomic hierarchy chart</p>
                </div>
            </div>
            <div class="figure-item">
                <h4>Diversity Analysis (PNG)</h4>
                <div class="figure-placeholder">
                    <i class="fas fa-chart-area" style="font-size: 3rem; color: var(--accent-green);"></i>
                    <p>Alpha & Beta diversity statistical plots</p>
                </div>
            </div>
            <div class="figure-item">
                <h4>Taxonomy Statistics (PDF)</h4>
                <div class="figure-placeholder">
                    <i class="fas fa-chart-bar" style="font-size: 3rem; color: var(--accent-cyan);"></i>
                    <p>Comprehensive taxonomic summary tables</p>
                </div>
            </div>
        </div>
        <div style="margin-top: 2rem; padding: 1rem; background: rgba(30, 106, 255, 0.1); border-radius: 8px; border: 1px solid rgba(30, 106, 255, 0.3);">
            <h4 style="color: var(--accent-blue); margin-bottom: 0.5rem;"><i class="fas fa-palette"></i> Custom Color Schemes Available:</h4>
            <div class="color-scheme-preview">
                <div class="color-dot" style="background: #1e6aff;"></div>
                <div class="color-dot" style="background: #00d4ff;"></div>
                <div class="color-dot" style="background: #00ff88;"></div>
                <div class="color-dot" style="background: #ff8c00;"></div>
                <div class="color-dot" style="background: #ff6b6b;"></div>
            </div>
            <p style="color: var(--text-gray); margin: 0.5rem 0 0 0; font-size: 0.9rem;">Ocean Blue • Cyan Depths • Marine Green • Deep Orange • Coral Red</p>
        </div>
    `;
    
    const figureModal = document.createElement('div');
    figureModal.className = 'modal';
    figureModal.style.display = 'block';
    figureModal.innerHTML = `
        <div class="modal-content" style="max-width: 900px;">
            <div class="modal-header">
                <h3>Publication Figures Gallery</h3>
                <span class="modal-close" onclick="this.closest('.modal').remove()">&times;</span>
            </div>
            <div class="modal-body">
                ${figureContent}
            </div>
        </div>
    `;
    
    document.body.appendChild(figureModal);
}

function customizeFigures() {
    showNotification('Opening figure customization panel...', 'info');
    
    const customizationContent = `
        <h3>Customize Publication Figures</h3>
        <div class="customization-panel">
            <div class="customization-options">
                <div class="customization-option">
                    <label for="colorScheme">Color Scheme</label>
                    <select id="colorScheme">
                        <option value="ocean">Ocean Blue Theme</option>
                        <option value="deep">Deep Sea Gradient</option>
                        <option value="marine">Marine Life Colors</option>
                        <option value="coral">Coral Reef Palette</option>
                        <option value="custom">Custom Colors</option>
                    </select>
                </div>
                <div class="customization-option">
                    <label for="resolution">Resolution (DPI)</label>
                    <select id="resolution">
                        <option value="150">150 DPI (Web)</option>
                        <option value="300" selected>300 DPI (Print)</option>
                        <option value="600">600 DPI (High Quality)</option>
                    </select>
                </div>
                <div class="customization-option">
                    <label for="format">Output Format</label>
                    <select id="format">
                        <option value="png">PNG (Raster)</option>
                        <option value="pdf" selected>PDF (Vector)</option>
                        <option value="svg">SVG (Vector)</option>
                        <option value="both">Both PNG & PDF</option>
                    </select>
                </div>
                <div class="customization-option">
                    <label for="fontSize">Font Size</label>
                    <input type="range" id="fontSize" min="8" max="16" value="12">
                    <span style="color: var(--text-gray); font-size: 0.8rem;">12pt</span>
                </div>
            </div>
            <div style="margin-top: 2rem; text-align: center;">
                <button class="btn btn-primary" onclick="generateCustomFigures()">
                    <i class="fas fa-magic"></i> Generate Custom Figures
                </button>
            </div>
        </div>
    `;
    
    const customModal = document.createElement('div');
    customModal.className = 'modal';
    customModal.style.display = 'block';
    customModal.innerHTML = `
        <div class="modal-content" style="max-width: 700px;">
            <div class="modal-header">
                <h3>Figure Customization</h3>
                <span class="modal-close" onclick="this.closest('.modal').remove()">&times;</span>
            </div>
            <div class="modal-body">
                ${customizationContent}
            </div>
        </div>
    `;
    
    document.body.appendChild(customModal);
    
    // Add font size slider interaction
    const fontSlider = customModal.querySelector('#fontSize');
    const fontDisplay = customModal.querySelector('span');
    fontSlider.addEventListener('input', (e) => {
        fontDisplay.textContent = e.target.value + 'pt';
    });
}

function generateCustomFigures() {
    showNotification('Generating custom figures with your settings...', 'info');
    
    // Simulate generation process
    let progress = 0;
    const progressInterval = setInterval(() => {
        progress += Math.random() * 20;
        if (progress >= 100) {
            progress = 100;
            clearInterval(progressInterval);
            showNotification('Custom figures generated successfully!', 'success');
            
            // Close the customization modal
            const modal = document.querySelector('.modal');
            if (modal) modal.remove();
        }
        console.log(`Figure generation progress: ${Math.round(progress)}%`);
    }, 500);
}

function downloadAllResults() {
    showNotification('Preparing complete analysis package...', 'info');
    
    // Show download progress
    let progress = 0;
    const progressInterval = setInterval(() => {
        progress += Math.random() * 15;
        if (progress >= 100) {
            progress = 100;
            clearInterval(progressInterval);
            showNotification('Complete analysis package downloaded successfully!', 'success');
        }
        console.log(`Download progress: ${Math.round(progress)}%`);
    }, 300);
}