// Main JavaScript functionality for DeepSea eDNA Explorer

document.addEventListener('DOMContentLoaded', function() {
    // Initialize all interactive features
    initNavigation();
    initStatCounters();
    initModal();
    initScrollAnimations();
    initUploadFeatures();
    initThemeToggle();
    initSpeciesCards();
});

// Navigation functionality
function initNavigation() {
    const navLinks = document.querySelectorAll('.nav-link');
    const sections = document.querySelectorAll('section[id]');

    // Smooth scrolling for navigation links
    navLinks.forEach(link => {
        link.addEventListener('click', function(e) {
            e.preventDefault();
            const targetId = this.getAttribute('href').substring(1);
            const targetSection = document.getElementById(targetId);
            
            if (targetSection) {
                targetSection.scrollIntoView({
                    behavior: 'smooth',
                    block: 'start'
                });
            }
        });
    });

    // Active link highlighting on scroll
    window.addEventListener('scroll', () => {
        let current = '';
        sections.forEach(section => {
            const sectionTop = section.offsetTop;
            const sectionHeight = section.clientHeight;
            if (pageYOffset >= sectionTop - 200) {
                current = section.getAttribute('id');
            }
        });

        navLinks.forEach(link => {
            link.classList.remove('active');
            if (link.getAttribute('href').substring(1) === current) {
                link.classList.add('active');
            }
        });
    });

    // Navbar background on scroll
    window.addEventListener('scroll', () => {
        const navbar = document.querySelector('.navbar');
        if (window.scrollY > 100) {
            navbar.style.background = 'rgba(10, 22, 40, 0.98)';
        } else {
            navbar.style.background = 'rgba(10, 22, 40, 0.95)';
        }
    });
}

// Animated statistics counters
function initStatCounters() {
    const statNumbers = document.querySelectorAll('.stat-number');
    
    const animateCounter = (element) => {
        const target = parseInt(element.getAttribute('data-target'));
        const duration = 2000; // 2 seconds
        const step = target / (duration / 16); // 60fps
        let current = 0;

        const timer = setInterval(() => {
            current += step;
            if (current >= target) {
                current = target;
                clearInterval(timer);
            }
            
            // Format numbers with commas for large numbers
            let displayValue = Math.floor(current);
            if (target > 1000) {
                displayValue = displayValue.toLocaleString();
            }
            
            // Add % for percentage values
            if (element.parentElement.querySelector('.stat-label').textContent.includes('Speed') ||
                element.parentElement.querySelector('.stat-label').textContent.includes('Improvement')) {
                displayValue += '%';
            } else if (target > 1000) {
                displayValue += '+';
            }
            
            element.textContent = displayValue;
        }, 16);
    };

    // Intersection Observer for counter animation
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

// Modal functionality
function initModal() {
    const modal = document.getElementById('uploadModal');
    const uploadBtn = document.getElementById('uploadBtn');
    const demoBtn = document.getElementById('demoBtn');
    const closeBtn = document.querySelector('.modal-close');

    // Open modal
    uploadBtn.addEventListener('click', () => {
        modal.style.display = 'block';
        document.body.style.overflow = 'hidden';
    });

    // Demo button functionality
    demoBtn.addEventListener('click', () => {
        // Redirect to analysis page or show demo
        window.location.href = 'analysis.html';
    });

    // Close modal
    closeBtn.addEventListener('click', () => {
        modal.style.display = 'none';
        document.body.style.overflow = 'auto';
    });

    // Close modal on outside click
    window.addEventListener('click', (e) => {
        if (e.target === modal) {
            modal.style.display = 'none';
            document.body.style.overflow = 'auto';
        }
    });

    // Close modal on escape key
    document.addEventListener('keydown', (e) => {
        if (e.key === 'Escape' && modal.style.display === 'block') {
            modal.style.display = 'none';
            document.body.style.overflow = 'auto';
        }
    });
}

// Scroll animations
function initScrollAnimations() {
    const animateElements = document.querySelectorAll('.tech-card, .stat-card, .problem-text, .performance-metrics');
    
    const observer = new IntersectionObserver((entries) => {
        entries.forEach(entry => {
            if (entry.isIntersecting) {
                entry.target.classList.add('animate-fade-in-up');
            }
        });
    }, {
        threshold: 0.1,
        rootMargin: '0px 0px -50px 0px'
    });

    animateElements.forEach(el => observer.observe(el));

    // Animate metric bars
    const metricBars = document.querySelectorAll('.metric-fill');
    const metricObserver = new IntersectionObserver((entries) => {
        entries.forEach(entry => {
            if (entry.isIntersecting && !entry.target.hasAttribute('data-animated')) {
                entry.target.setAttribute('data-animated', 'true');
                const width = entry.target.style.width;
                entry.target.style.width = '0%';
                setTimeout(() => {
                    entry.target.style.width = width;
                }, 100);
            }
        });
    });

    metricBars.forEach(bar => metricObserver.observe(bar));
}

// Upload functionality
function initUploadFeatures() {
    const uploadArea = document.getElementById('uploadArea');
    const fileInput = document.getElementById('fileInput');

    // Drag and drop functionality
    uploadArea.addEventListener('dragover', (e) => {
        e.preventDefault();
        uploadArea.style.borderColor = 'var(--accent-blue)';
        uploadArea.style.background = 'rgba(30, 106, 255, 0.1)';
    });

    uploadArea.addEventListener('dragleave', (e) => {
        e.preventDefault();
        uploadArea.style.borderColor = 'rgba(30, 106, 255, 0.3)';
        uploadArea.style.background = 'transparent';
    });

    uploadArea.addEventListener('drop', (e) => {
        e.preventDefault();
        uploadArea.style.borderColor = 'rgba(30, 106, 255, 0.3)';
        uploadArea.style.background = 'transparent';
        
        const files = e.dataTransfer.files;
        if (files.length > 0) {
            handleFileUpload(files[0]);
        }
    });

    // Click to upload
    uploadArea.addEventListener('click', () => {
        fileInput.click();
    });

    // File input change
    fileInput.addEventListener('change', (e) => {
        if (e.target.files.length > 0) {
            handleFileUpload(e.target.files[0]);
        }
    });
}

// Handle file upload
function handleFileUpload(file) {
    // Validate file type
    const allowedTypes = ['.fasta', '.fastq', '.fa', '.fq'];
    const fileExtension = '.' + file.name.split('.').pop().toLowerCase();
    
    if (!allowedTypes.includes(fileExtension)) {
        showNotification('Please upload a valid eDNA sequence file (.fasta, .fastq, .fa, .fq)', 'error');
        return;
    }

    // Validate file size (max 100MB)
    if (file.size > 100 * 1024 * 1024) {
        showNotification('File size too large. Please upload a file smaller than 100MB.', 'error');
        return;
    }

    // Show upload progress
    showUploadProgress(file);
    
    // Simulate file processing
    setTimeout(() => {
        hideUploadProgress();
        showNotification('File uploaded successfully! Redirecting to analysis...', 'success');
        setTimeout(() => {
            window.location.href = 'analysis.html';
        }, 2000);
    }, 3000);
}

// Show upload progress
function showUploadProgress(file) {
    const uploadArea = document.getElementById('uploadArea');
    uploadArea.innerHTML = `
        <div class="upload-progress">
            <i class="fas fa-dna loading"></i>
            <h3>Processing ${file.name}</h3>
            <div class="progress-bar">
                <div class="progress-fill"></div>
            </div>
            <p>Preparing your eDNA sample for analysis...</p>
        </div>
    `;

    // Animate progress bar
    const progressFill = uploadArea.querySelector('.progress-fill');
    let progress = 0;
    const progressTimer = setInterval(() => {
        progress += Math.random() * 15;
        if (progress >= 100) {
            progress = 100;
            clearInterval(progressTimer);
        }
        progressFill.style.width = progress + '%';
    }, 200);
}

// Hide upload progress
function hideUploadProgress() {
    const uploadArea = document.getElementById('uploadArea');
    uploadArea.innerHTML = `
        <i class="fas fa-cloud-upload-alt"></i>
        <p>Drag and drop your eDNA sequence file here or click to browse</p>
        <input type="file" id="fileInput" accept=".fasta,.fastq,.fa,.fq" hidden>
        <button class="btn btn-outline" onclick="document.getElementById('fileInput').click()">
            Choose File
        </button>
    `;
}

// Notification system
function showNotification(message, type = 'info') {
    const notification = document.createElement('div');
    notification.className = `notification notification-${type}`;
    notification.innerHTML = `
        <i class="fas fa-${type === 'success' ? 'check-circle' : type === 'error' ? 'exclamation-circle' : 'info-circle'}"></i>
        <span>${message}</span>
        <button class="notification-close">&times;</button>
    `;

    document.body.appendChild(notification);

    // Add styles dynamically
    const style = document.createElement('style');
    style.textContent = `
        .notification {
            position: fixed;
            top: 100px;
            right: 20px;
            background: var(--secondary-dark);
            border: 1px solid ${type === 'success' ? 'var(--accent-green)' : type === 'error' ? 'var(--accent-orange)' : 'var(--accent-blue)'};
            border-radius: var(--border-radius);
            padding: 1rem 1.5rem;
            display: flex;
            align-items: center;
            gap: 0.75rem;
            z-index: 3000;
            animation: slideInRight 0.3s ease-out;
            max-width: 400px;
            box-shadow: var(--shadow-card);
        }
        
        .notification i {
            color: ${type === 'success' ? 'var(--accent-green)' : type === 'error' ? 'var(--accent-orange)' : 'var(--accent-blue)'};
            font-size: 1.2rem;
        }
        
        .notification-close {
            background: none;
            border: none;
            color: var(--text-gray);
            font-size: 1.2rem;
            cursor: pointer;
            padding: 0;
            margin-left: auto;
        }
        
        @keyframes slideInRight {
            from {
                transform: translateX(100%);
                opacity: 0;
            }
            to {
                transform: translateX(0);
                opacity: 1;
            }
        }
        
        .upload-progress {
            text-align: center;
            padding: 2rem;
        }
        
        .upload-progress h3 {
            margin: 1rem 0;
            color: var(--text-white);
        }
        
        .progress-bar {
            width: 100%;
            height: 8px;
            background: rgba(255, 255, 255, 0.1);
            border-radius: 4px;
            margin: 1rem 0;
            overflow: hidden;
        }
        
        .progress-fill {
            height: 100%;
            background: var(--gradient-primary);
            width: 0%;
            transition: width 0.3s ease;
            border-radius: 4px;
        }
    `;
    
    if (!document.head.querySelector('style[data-notifications]')) {
        style.setAttribute('data-notifications', 'true');
        document.head.appendChild(style);
    }

    // Close notification
    notification.querySelector('.notification-close').addEventListener('click', () => {
        notification.remove();
    });

    // Auto remove after 5 seconds
    setTimeout(() => {
        if (notification.parentNode) {
            notification.remove();
        }
    }, 5000);
}

// Theme toggle functionality
function initThemeToggle() {
    const themeToggle = document.querySelector('.nav-toggle');
    let isDark = true;

    themeToggle.addEventListener('click', () => {
        isDark = !isDark;
        
        if (isDark) {
            themeToggle.innerHTML = '<i class="fas fa-moon"></i>';
            document.documentElement.style.setProperty('--primary-dark', '#0a1628');
            document.documentElement.style.setProperty('--secondary-dark', '#1e2936');
        } else {
            themeToggle.innerHTML = '<i class="fas fa-sun"></i>';
            document.documentElement.style.setProperty('--primary-dark', '#f8fafc');
            document.documentElement.style.setProperty('--secondary-dark', '#e2e8f0');
        }
    });
}

// Utility functions
function debounce(func, wait) {
    let timeout;
    return function executedFunction(...args) {
        const later = () => {
            clearTimeout(timeout);
            func(...args);
        };
        clearTimeout(timeout);
        timeout = setTimeout(later, wait);
    };
}

// Performance optimization
const debouncedScroll = debounce(() => {
    // Scroll-based animations
}, 100);

window.addEventListener('scroll', debouncedScroll);

// Species data
const speciesData = {
    'vampire-squid': {
        name: 'Vampire Squid',
        scientificName: 'Vampyroteuthis infernalis',
        description: 'A unique cephalopod that inhabits the deep ocean with light-producing organs and webbed arms. Despite its name, it doesn\'t suck blood but feeds on marine snow and detritus that drifts down from the surface.',
        depth: '500-3000m',
        habitat: 'Oxygen Minimum Zones of Pacific Ocean',
        size: 'Up to 30cm',
        conservationStatus: 'Least Concern',
        diet: 'Marine snow, detritus',
        lifespan: '8+ years',
        specialFeatures: 'Bioluminescent photophores, inside-out defense mechanism',
        habitatDescription: 'Lives in the oxygen minimum zone where oxygen levels are too low for most other organisms to survive. This unique adaptation allows it to thrive in one of the ocean\'s most challenging environments.'
    },
    'anglerfish': {
        name: 'Anglerfish',
        scientificName: 'Lophiiformes',
        description: 'Known for its bioluminescent lure that attracts prey in the dark depths. Female anglerfish are much larger than males, with some species showing extreme sexual dimorphism.',
        depth: '200-2000m',
        habitat: 'Deep ocean waters worldwide',
        size: 'Up to 1.2m (females)',
        conservationStatus: 'Various species status',
        diet: 'Small fish, crustaceans',
        lifespan: '25+ years',
        specialFeatures: 'Bioluminescent lure, extreme sexual dimorphism',
        habitatDescription: 'Found in deep waters across all oceans, using their bioluminescent lure to attract prey in the pitch-black depths where sunlight cannot penetrate.'
    },
    'giant-tube-worm': {
        name: 'Giant Tube Worm',
        scientificName: 'Riftia pachyptila',
        description: 'Thrives near hydrothermal vents and relies on symbiotic bacteria for nutrition. These remarkable creatures have no mouth or stomach, instead depending entirely on chemosynthetic bacteria.',
        depth: '2000-2800m',
        habitat: 'Hydrothermal vents of Pacific Ocean',
        size: 'Up to 2.4m',
        conservationStatus: 'Not Evaluated',
        diet: 'Chemosynthetic bacteria (symbiotic)',
        lifespan: '300+ years',
        specialFeatures: 'No digestive system, rapid growth, heat tolerance',
        habitatDescription: 'Lives exclusively around hydrothermal vents where they form dense colonies. These vents provide the chemical energy needed by their symbiotic bacteria to produce nutrients.'
    },
    'deep-sea-jellyfish': {
        name: 'Deep Sea Jellyfish',
        scientificName: 'Atolla jellyfish',
        description: 'Transparent gelatinous creature with bioluminescent capabilities for defense. When threatened, it creates a "burglar alarm" - a ring of light to attract predators to whatever is attacking it.',
        depth: '1000-4000m',
        habitat: 'Deep ocean waters worldwide',
        size: '10-20cm diameter',
        conservationStatus: 'Not Evaluated',
        diet: 'Small fish, plankton, other jellyfish',
        lifespan: '1-2 years',
        specialFeatures: 'Bioluminescent alarm system, transparent body',
        habitatDescription: 'Inhabits the midnight zone of oceans worldwide, using its unique bioluminescent defense mechanism to survive in the predator-rich deep sea environment.'
    },
    'dumbo-octopus': {
        name: 'Dumbo Octopus',
        scientificName: 'Grimpoteuthis',
        description: 'Deepest living octopus with ear-like fins that propel it through the water. Named after Disney\'s Dumbo elephant due to its distinctive ear-like fins.',
        depth: '3000-7000m',
        habitat: 'Abyssal plains worldwide',
        size: '20-30cm',
        conservationStatus: 'Not Evaluated',
        diet: 'Worms, crustaceans, copepods',
        lifespan: '3-5 years',
        specialFeatures: 'Ear-like fins, deepest-living octopus, gelatinous body',
        habitatDescription: 'Lives on or near the ocean floor in the abyssal zone, using its ear-like fins to gracefully navigate the extreme depths where pressure can exceed 600 times that at sea level.'
    },
    'barreleye-fish': {
        name: 'Barreleye Fish',
        scientificName: 'Macropinna microstoma',
        description: 'Features a transparent head and tubular eyes that can rotate upward. This unique adaptation allows it to see through its own transparent head to spot prey silhouetted against the faint light above.',
        depth: '600-800m',
        habitat: 'North Pacific Ocean',
        size: '15cm',
        conservationStatus: 'Not Evaluated',
        diet: 'Small fish, jellyfish tentacles',
        lifespan: '5-10 years',
        specialFeatures: 'Transparent head, rotating tubular eyes, light-sensing ability',
        habitatDescription: 'Lives in the twilight zone where some sunlight still penetrates. Its transparent head and rotating eyes are perfectly adapted to spot prey and predators in this dimly lit environment.'
    }
};

// Initialize species functionality
function initSpeciesCards() {
    const speciesCards = document.querySelectorAll('.species-card');
    const speciesModal = document.getElementById('speciesModal');
    const speciesModalClose = document.getElementById('speciesModalClose');
    
    speciesCards.forEach(card => {
        card.addEventListener('click', () => {
            const speciesId = card.getAttribute('data-species');
            showSpeciesModal(speciesId);
        });
    });
    
    // Close modal functionality
    if (speciesModalClose) {
        speciesModalClose.addEventListener('click', () => {
            speciesModal.style.display = 'none';
            document.body.style.overflow = 'auto';
        });
    }
    
    // Close modal on outside click
    window.addEventListener('click', (e) => {
        if (e.target === speciesModal) {
            speciesModal.style.display = 'none';
            document.body.style.overflow = 'auto';
        }
    });
    
    // Close modal on escape key
    document.addEventListener('keydown', (e) => {
        if (e.key === 'Escape' && speciesModal.style.display === 'block') {
            speciesModal.style.display = 'none';
            document.body.style.overflow = 'auto';
        }
    });
}

// Show species modal with details
function showSpeciesModal(speciesId) {
    const species = speciesData[speciesId];
    if (!species) return;
    
    const modal = document.getElementById('speciesModal');
    const modalTitle = document.getElementById('speciesModalTitle');
    const modalContent = document.getElementById('speciesDetailContent');
    
    modalTitle.textContent = species.name;
    
    modalContent.innerHTML = `
        <div class="species-detail-header">
            <div class="species-detail-image">
                <div class="species-placeholder ${speciesId}" style="height: 100%; border-radius: 8px;"></div>
            </div>
            <div class="species-detail-info">
                <h2>${species.name}</h2>
                <p class="species-scientific-name">${species.scientificName}</p>
                <p class="species-detail-description">${species.description}</p>
            </div>
        </div>
        
        <div class="species-facts">
            <div class="fact-item">
                <div class="fact-label">Depth Range</div>
                <div class="fact-value">${species.depth}</div>
            </div>
            <div class="fact-item">
                <div class="fact-label">Size</div>
                <div class="fact-value">${species.size}</div>
            </div>
            <div class="fact-item">
                <div class="fact-label">Diet</div>
                <div class="fact-value">${species.diet}</div>
            </div>
            <div class="fact-item">
                <div class="fact-label">Lifespan</div>
                <div class="fact-value">${species.lifespan}</div>
            </div>
            <div class="fact-item">
                <div class="fact-label">Conservation Status</div>
                <div class="fact-value">${species.conservationStatus}</div>
            </div>
            <div class="fact-item">
                <div class="fact-label">Special Features</div>
                <div class="fact-value">${species.specialFeatures}</div>
            </div>
        </div>
        
        <div class="species-habitat">
            <h4><i class="fas fa-map-marker-alt"></i> Habitat</h4>
            <div class="habitat-description">
                <strong>${species.habitat}</strong><br>
                ${species.habitatDescription}
            </div>
        </div>
        
        <button class="learn-more-btn" onclick="window.open('https://en.wikipedia.org/wiki/${encodeURIComponent(species.scientificName)}', '_blank')">
            <i class="fas fa-external-link-alt"></i>
            Learn More
        </button>
    `;
    
    modal.style.display = 'block';
    document.body.style.overflow = 'hidden';
}

// Preload critical resources
function preloadResources() {
    const criticalImages = [
        // Add any critical images here
    ];
    
    criticalImages.forEach(src => {
        const link = document.createElement('link');
        link.rel = 'preload';
        link.as = 'image';
        link.href = src;
        document.head.appendChild(link);
    });
}

// Initialize on page load
preloadResources();