/**
 * PharmaSight™ - Loading Screens & Educational Slides
 * Manages loading animations and displays educational content during wait times
 */

class LoadingScreen {
    constructor() {
        this.loadingElement = document.getElementById('loading-screen');
        this.slideElement = document.querySelector('.educational-slide .slide-text');
        this.currentSlideIndex = 0;
        this.slideInterval = null;

        // Educational facts about drug discovery
        this.educationalSlides = [
            "Did you know? Drug discovery typically takes 10-15 years and costs $2.6 billion",
            "Only 1 in 5,000 compounds makes it from discovery to FDA approval",
            "AI can reduce drug discovery time by up to 75% through computational screening",
            "The human body has over 20,000 proteins that could be potential drug targets",
            "Quantum chemistry calculations can predict molecular properties with 95%+ accuracy",
            "PharmaSight analyzes millions of molecular combinations in minutes",
            "ADMET prediction helps eliminate unsafe compounds early in development",
            "Protein-ligand docking simulations predict drug efficacy before synthesis",
            "The blood-brain barrier blocks 98% of small molecules from entering the brain",
            "Machine learning models can predict drug toxicity from molecular structure",
            "Retrosynthesis AI suggests optimal routes to synthesize new compounds",
            "Generative chemistry creates novel molecules never seen in nature",
            "Over 50% of approved drugs work by binding to G-protein coupled receptors",
            "Drug metabolism prediction prevents expensive clinical trial failures",
            "Computational chemistry saves millions in lab costs per compound",
            "Molecular fingerprints help find similar compounds with better properties",
            "Quantum tunneling affects how drugs interact at the atomic level",
            "AI can identify drug repurposing opportunities in existing medications",
            "Bioavailability determines how much of a drug reaches its target",
            "PharmaSight integrates 15+ computational chemistry tools in one platform"
        ];
    }

    show() {
        if (this.loadingElement) {
            this.loadingElement.classList.remove('hidden');
            this.startSlideShow();
        }
    }

    hide(delay = 500) {
        setTimeout(() => {
            if (this.loadingElement) {
                this.loadingElement.classList.add('hidden');
                this.stopSlideShow();
            }
        }, delay);
    }

    startSlideShow() {
        // Show first slide immediately
        this.updateSlide();

        // Rotate slides every 4 seconds
        this.slideInterval = setInterval(() => {
            this.updateSlide();
        }, 4000);
    }

    stopSlideShow() {
        if (this.slideInterval) {
            clearInterval(this.slideInterval);
            this.slideInterval = null;
        }
    }

    updateSlide() {
        if (!this.slideElement) return;

        // Fade out
        this.slideElement.style.opacity = '0';
        this.slideElement.style.transform = 'translateY(10px)';

        setTimeout(() => {
            // Update text
            this.slideElement.textContent = this.educationalSlides[this.currentSlideIndex];
            this.currentSlideIndex = (this.currentSlideIndex + 1) % this.educationalSlides.length;

            // Fade in
            this.slideElement.style.transition = 'all 0.5s ease';
            this.slideElement.style.opacity = '1';
            this.slideElement.style.transform = 'translateY(0)';
        }, 300);
    }

    updateText(text) {
        const loadingText = document.querySelector('.loading-text');
        if (loadingText) {
            loadingText.textContent = text;
        }
    }

    showCustomSlide(text) {
        if (this.slideElement) {
            this.slideElement.style.opacity = '0';
            setTimeout(() => {
                this.slideElement.textContent = text;
                this.slideElement.style.opacity = '1';
            }, 300);
        }
    }
}

// Simulated API Loading Screen
class APILoadingManager {
    constructor() {
        this.loadingScreen = new LoadingScreen();
        this.processingStages = [
            { text: 'Analyzing molecular structure...', duration: 1000 },
            { text: 'Calculating ADMET properties...', duration: 1200 },
            { text: 'Running quantum chemistry calculations...', duration: 1500 },
            { text: 'Predicting toxicity profile...', duration: 1000 },
            { text: 'Generating synthesis routes...', duration: 1300 },
            { text: 'Finalizing results...', duration: 800 }
        ];
    }

    async simulateAPICall(endpoint, data) {
        this.loadingScreen.show();

        // Simulate processing stages
        for (const stage of this.processingStages) {
            this.loadingScreen.updateText(stage.text);
            await this.delay(stage.duration);
        }

        // Make actual API call
        try {
            const response = await fetch(endpoint, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(data)
            });

            const result = await response.json();
            this.loadingScreen.hide();
            return result;
        } catch (error) {
            this.loadingScreen.updateText('Error processing request');
            this.loadingScreen.hide(2000);
            throw error;
        }
    }

    delay(ms) {
        return new Promise(resolve => setTimeout(resolve, ms));
    }
}

// Initialize loading screen manager
const loadingManager = new APILoadingManager();

// Hide loading screen when page loads
window.addEventListener('load', () => {
    const loadingScreen = new LoadingScreen();
    loadingScreen.hide(1000);
});

// Export for use in other modules
if (typeof module !== 'undefined' && module.exports) {
    module.exports = { LoadingScreen, APILoadingManager };
}
