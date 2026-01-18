/**
 * PharmaSight™ - Main Application Logic
 * Handles scroll animations, navigation, and general interactivity
 */

class PharmaSightApp {
    constructor() {
        this.init();
        this.setupScrollAnimations();
        this.setupNavigation();
        this.setupInteractions();
    }

    init() {
        console.log('PharmaSight™ Platform Initialized');

        // Add smooth scroll behavior
        document.querySelectorAll('a[href^="#"]').forEach(anchor => {
            anchor.addEventListener('click', (e) => {
                e.preventDefault();
                const target = document.querySelector(anchor.getAttribute('href'));
                if (target) {
                    target.scrollIntoView({ behavior: 'smooth', block: 'start' });
                }
            });
        });
    }

    setupScrollAnimations() {
        // Intersection Observer for fade-up animations
        const observerOptions = {
            root: null,
            rootMargin: '0px',
            threshold: 0.1
        };

        const observer = new IntersectionObserver((entries) => {
            entries.forEach(entry => {
                if (entry.isIntersecting) {
                    entry.target.classList.add('animated');
                }
            });
        }, observerOptions);

        // Observe all elements with data-animation attribute
        document.querySelectorAll('[data-animation]').forEach(el => {
            observer.observe(el);
        });

        // Parallax effect on scroll
        let ticking = false;
        window.addEventListener('scroll', () => {
            if (!ticking) {
                window.requestAnimationFrame(() => {
                    this.handleParallax();
                    ticking = false;
                });
                ticking = true;
            }
        });
    }

    handleParallax() {
        const scrolled = window.pageYOffset;
        const parallaxElements = document.querySelectorAll('.parallax-layer');

        parallaxElements.forEach(el => {
            const speed = el.dataset.speed || 0.5;
            const yPos = -(scrolled * speed);
            el.style.transform = `translateY(${yPos}px)`;
        });

        // Update navigation background on scroll
        const nav = document.querySelector('.nav-glass');
        if (nav) {
            if (scrolled > 50) {
                nav.style.background = 'rgba(10, 10, 15, 0.95)';
                nav.style.boxShadow = '0 4px 20px rgba(0, 0, 0, 0.3)';
            } else {
                nav.style.background = 'rgba(10, 10, 15, 0.8)';
                nav.style.boxShadow = 'none';
            }
        }
    }

    setupNavigation() {
        // Active section highlighting
        const sections = document.querySelectorAll('section[id]');
        const navLinks = document.querySelectorAll('.nav-link');

        const highlightNav = () => {
            const scrollY = window.pageYOffset;

            sections.forEach(section => {
                const sectionHeight = section.offsetHeight;
                const sectionTop = section.offsetTop - 100;
                const sectionId = section.getAttribute('id');

                if (scrollY > sectionTop && scrollY <= sectionTop + sectionHeight) {
                    navLinks.forEach(link => {
                        link.classList.remove('active');
                        if (link.getAttribute('href') === `#${sectionId}`) {
                            link.classList.add('active');
                        }
                    });
                }
            });
        };

        window.addEventListener('scroll', highlightNav);
    }

    setupInteractions() {
        // Button ripple effect
        document.querySelectorAll('.btn-primary, .btn-secondary').forEach(button => {
            button.addEventListener('click', function (e) {
                const ripple = document.createElement('span');
                const rect = this.getBoundingClientRect();
                const size = Math.max(rect.width, rect.height);
                const x = e.clientX - rect.left - size / 2;
                const y = e.clientY - rect.top - size / 2;

                ripple.style.cssText = `
                    position: absolute;
                    width: ${size}px;
                    height: ${size}px;
                    border-radius: 50%;
                    background: rgba(255, 255, 255, 0.4);
                    left: ${x}px;
                    top: ${y}px;
                    pointer-events: none;
                    animation: ripple-effect 0.6s ease-out;
                `;

                this.appendChild(ripple);

                setTimeout(() => ripple.remove(), 600);
            });
        });

        // Card hover 3D effect
        document.querySelectorAll('.feature-card, .glass-card').forEach(card => {
            card.addEventListener('mousemove', (e) => {
                const rect = card.getBoundingClientRect();
                const x = e.clientX - rect.left;
                const y = e.clientY - rect.top;

                const centerX = rect.width / 2;
                const centerY = rect.height / 2;

                const rotateX = (y - centerY) / 20;
                const rotateY = (centerX - x) / 20;

                card.style.transform = `perspective(1000px) rotateX(${rotateX}deg) rotateY(${rotateY}deg) translateY(-4px)`;
            });

            card.addEventListener('mouseleave', () => {
                card.style.transform = 'perspective(1000px) rotateX(0) rotateY(0)';
            });
        });

        // Counter animation for stats
        this.animateCounters();
    }

    animateCounters() {
        const counters = document.querySelectorAll('.stat-number');
        const duration = 2000; // 2 seconds

        const observerOptions = {
            root: null,
            threshold: 0.5
        };

        const observer = new IntersectionObserver((entries) => {
            entries.forEach(entry => {
                if (entry.isIntersecting && !entry.target.dataset.counted) {
                    entry.target.dataset.counted = 'true';
                    const target = entry.target;
                    const text = target.textContent;

                    // Extract number and suffix
                    const match = text.match(/([0-9,.]+)(.*)$/);
                    if (!match) return;

                    const numStr = match[1].replace(/,/g, '');
                    const suffix = match[2];
                    const isDecimal = numStr.includes('.');

                    if (isDecimal) {
                        // Handle decimal numbers (like 98.5%)
                        const finalNum = parseFloat(numStr);
                        let current = 0;
                        const increment = finalNum / (duration / 16);

                        const timer = setInterval(() => {
                            current += increment;
                            if (current >= finalNum) {
                                current = finalNum;
                                clearInterval(timer);
                            }
                            target.textContent = current.toFixed(1) + suffix;
                        }, 16);
                    } else {
                        // Handle integers (like 10,000+)
                        const finalNum = parseInt(numStr);
                        let current = 0;
                        const increment = finalNum / (duration / 16);

                        const timer = setInterval(() => {
                            current += increment;
                            if (current >= finalNum) {
                                current = finalNum;
                                clearInterval(timer);
                            }
                            target.textContent = Math.floor(current).toLocaleString() + suffix;
                        }, 16);
                    }
                }
            });
        }, observerOptions);

        counters.forEach(counter => observer.observe(counter));
    }

    // API Integration Methods
    async fetchWithLoading(endpoint, options = {}) {
        const loadingManager = new APILoadingManager();
        try {
            return await loadingManager.simulateAPICall(endpoint, options.body ? JSON.parse(options.body) : {});
        } catch (error) {
            console.error('API Error:', error);
            this.showNotification('An error occurred. Please try again.', 'error');
        }
    }

    showNotification(message, type = 'info') {
        const notification = document.createElement('div');
        notification.className = `notification notification-${type}`;
        notification.style.cssText = `
            position: fixed;
            top: 20px;
            right: 20px;
            padding: 1rem 2rem;
            background: ${type === 'error' ? 'rgba(255, 0, 110, 0.9)' : 'rgba(0, 212, 255, 0.9)'};
            color: white;
            border-radius: 1rem;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.3);
            backdrop-filter: blur(10px);
            z-index: 10000;
            animation: slideInRight 0.3s ease;
        `;
        notification.textContent = message;

        document.body.appendChild(notification);

        setTimeout(() => {
            notification.style.animation = 'slideOutRight 0.3s ease';
            setTimeout(() => notification.remove(), 300);
        }, 3000);
    }
}

// Add CSS animations for notifications
const style = document.createElement('style');
style.textContent = `
    @keyframes slideInRight {
        from {
            transform: translateX(400px);
            opacity: 0;
        }
        to {
            transform: translateX(0);
            opacity: 1;
        }
    }

    @keyframes slideOutRight {
        from {
            transform: translateX(0);
            opacity: 1;
        }
        to {
            transform: translateX(400px);
            opacity: 0;
        }
    }

    @keyframes ripple-effect {
        0% {
            transform: scale(0);
            opacity: 1;
        }
        100% {
            transform: scale(2);
            opacity: 0;
        }
    }
`;
document.head.appendChild(style);

// Initialize app when DOM is ready
document.addEventListener('DOMContentLoaded', () => {
    window.pharmasightApp = new PharmaSightApp();
});

// Expose API for external use
window.PharmaSight = {
    showNotification: (message, type) => {
        if (window.pharmasightApp) {
            window.pharmasightApp.showNotification(message, type);
        }
    },
    fetchWithLoading: (endpoint, options) => {
        if (window.pharmasightApp) {
            return window.pharmasightApp.fetchWithLoading(endpoint, options);
        }
    }
};
