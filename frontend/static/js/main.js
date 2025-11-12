// PharmaSight™ Frontend JavaScript

// API Configuration
// Support both local development and production deployment
const API_BASE = window.location.origin.includes('localhost') 
    ? 'http://localhost:8080' 
    : window.location.origin;

// Research Engine API - direct access for now
// TODO: Route through API Gateway once configured
const RESEARCH_API = window.location.origin.includes('localhost')
    ? 'http://localhost:8006'
    : `${window.location.protocol}//${window.location.hostname}:8006`;

// Utility Functions
function scrollToSection(sectionId) {
    document.getElementById(sectionId).scrollIntoView({ 
        behavior: 'smooth' 
    });
}

// Navigation
document.addEventListener('DOMContentLoaded', () => {
    // Active nav link on scroll
    const sections = document.querySelectorAll('section[id]');
    const navLinks = document.querySelectorAll('.nav-link');
    
    window.addEventListener('scroll', () => {
        let current = '';
        sections.forEach(section => {
            const sectionTop = section.offsetTop;
            const sectionHeight = section.clientHeight;
            if (window.pageYOffset >= sectionTop - 200) {
                current = section.getAttribute('id');
            }
        });
        
        navLinks.forEach(link => {
            link.classList.remove('active');
            if (link.getAttribute('href') === `#${current}`) {
                link.classList.add('active');
            }
        });
    });
    
    // Load initial data
    loadArticles();
    loadAnalogs();
    loadSystemHealth();
    initializeCharts();
});

// Research Articles
async function loadArticles() {
    try {
        const response = await fetch(`${RESEARCH_API}/articles/all`);
        const data = await response.json();
        
        if (data.success) {
            displayArticles(data.articles);
            populateYearFilter(data.articles);
            updateArticleCount(data.count);
        }
    } catch (error) {
        console.error('Error loading articles:', error);
        document.getElementById('articlesGrid').innerHTML = `
            <div class="error-message">
                <i class="fas fa-exclamation-triangle"></i>
                <p>Unable to load articles. Please check if the research engine is running.</p>
            </div>
        `;
    }
}

function displayArticles(articles) {
    const grid = document.getElementById('articlesGrid');
    
    if (articles.length === 0) {
        grid.innerHTML = '<p class="no-results">No articles found</p>';
        return;
    }
    
    grid.innerHTML = articles.map(article => `
        <div class="article-card">
            <div class="article-header">
                <h3>${article.title}</h3>
                <span class="article-year">${article.year || 'N/A'}</span>
            </div>
            <p class="article-authors">${article.authors ? article.authors.join(', ') : 'Unknown'}</p>
            <p class="article-journal"><i class="fas fa-book"></i> ${article.journal || 'Unknown'}</p>
            <div class="article-meta">
                ${article.doi ? `<span class="meta-tag"><i class="fas fa-link"></i> DOI: ${article.doi}</span>` : ''}
                ${article.pmid ? `<span class="meta-tag"><i class="fas fa-hashtag"></i> PMID: ${article.pmid}</span>` : ''}
            </div>
            ${article.doi ? `
                <a href="https://doi.org/${article.doi}" target="_blank" class="btn btn-outline btn-small">
                    <i class="fas fa-external-link-alt"></i> View Article
                </a>
            ` : ''}
        </div>
    `).join('');
}

function populateYearFilter(articles) {
    const years = [...new Set(articles.map(a => a.year).filter(Boolean))].sort((a, b) => b - a);
    const select = document.getElementById('filterYear');
    
    years.forEach(year => {
        const option = document.createElement('option');
        option.value = year;
        option.textContent = year;
        select.appendChild(option);
    });
}

function updateArticleCount(count) {
    const counter = document.getElementById('articleCount');
    if (counter) {
        animateCounter(counter, count);
    }
}

// Analog Discoveries
async function loadAnalogs() {
    try {
        const response = await fetch(`${RESEARCH_API}/discoveries/all`);
        const data = await response.json();
        
        if (data.success) {
            displayAnalogStats(data);
            updateAnalogCount(data.total_analogs);
        }
    } catch (error) {
        console.error('Error loading analogs:', error);
    }
}

function displayAnalogStats(data) {
    // Update stats in the analogs section
    console.log('Analog discoveries loaded:', data);
}

function updateAnalogCount(count) {
    const counter = document.getElementById('analogCount');
    if (counter) {
        animateCounter(counter, count);
    }
}

// System Health
async function loadSystemHealth() {
    try {
        const response = await fetch(`${API_BASE}/health`);
        const data = await response.json();
        
        displaySystemHealth(data);
    } catch (error) {
        console.error('Error loading system health:', error);
    }
}

function displaySystemHealth(health) {
    const container = document.getElementById('systemHealth');
    if (!container) return;
    
    const services = health.services || {};
    
    container.innerHTML = Object.entries(services).map(([name, status]) => `
        <div class="health-item">
            <span>${name.replace(/-/g, ' ').toUpperCase()}</span>
            <span class="health-status ${status.status === 'healthy' ? 'healthy' : 'unhealthy'}">
                ${status.status}
            </span>
        </div>
    `).join('');
}

// Research Cycle
async function runResearchCycle() {
    const btn = event.target;
    btn.disabled = true;
    btn.innerHTML = '<i class="fas fa-spinner fa-spin"></i> Running...';
    
    try {
        const response = await fetch(`${RESEARCH_API}/research/run-cycle`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                goals: [
                    "psilocybin 5-HT2A receptor depression",
                    "ketamine NMDA antagonist",
                    "MDMA PTSD therapy",
                    "kava lactones anxiolytic",
                    "muscimol GABA-A agonist"
                ]
            })
        });
        
        const data = await response.json();
        
        if (data.success) {
            alert(`Research cycle complete!\nArticles found: ${data.summary.articles_found || 0}`);
            loadArticles(); // Reload articles
        }
    } catch (error) {
        console.error('Error running research cycle:', error);
        alert('Error running research cycle. Please try again.');
    } finally {
        btn.disabled = false;
        btn.innerHTML = '<i class="fas fa-sync"></i> Run Research Cycle';
    }
}

// Search and Filter
document.getElementById('searchArticles')?.addEventListener('input', (e) => {
    // Implement search functionality
    console.log('Search:', e.target.value);
});

document.getElementById('filterYear')?.addEventListener('change', (e) => {
    // Implement filter functionality
    console.log('Filter year:', e.target.value);
});

// Feature Actions
function openCompoundAnalysis() {
    alert('Compound Analysis feature - Connect to compound-service API');
}

function openAnalogGeneration() {
    alert('Analog Generation feature - Connect to analog-service API');
}

function openResearchEngine() {
    scrollToSection('research');
}

function openPKPDSimulation() {
    alert('PKPD Simulation feature - Connect to ml-service API');
}

function viewPublicRegistry() {
    window.open('/PUBLIC_ANALOG_DISCOVERY_REGISTRY.md', '_blank');
}

// Counter Animation
function animateCounter(element, target) {
    const duration = 2000;
    const start = 0;
    const increment = target / (duration / 16);
    let current = start;
    
    const timer = setInterval(() => {
        current += increment;
        if (current >= target) {
            element.textContent = target;
            clearInterval(timer);
        } else {
            element.textContent = Math.floor(current);
        }
    }, 16);
}

// Charts
function initializeCharts() {
    // Database Statistics Chart
    const dbCtx = document.getElementById('databaseChart');
    if (dbCtx) {
        new Chart(dbCtx, {
            type: 'bar',
            data: {
                labels: ['Compounds', 'Analogs', 'Articles', 'Receptors'],
                datasets: [{
                    label: 'Count',
                    data: [500, 119, 26, 70],
                    backgroundColor: [
                        'rgba(0, 102, 204, 0.8)',
                        'rgba(0, 168, 232, 0.8)',
                        'rgba(0, 217, 255, 0.8)',
                        'rgba(0, 230, 118, 0.8)'
                    ]
                }]
            },
            options: {
                responsive: true,
                plugins: {
                    legend: { display: false }
                },
                scales: {
                    y: { beginAtZero: true }
                }
            }
        });
    }
    
    // Discovery Distribution Chart
    const discCtx = document.getElementById('discoveryChart');
    if (discCtx) {
        new Chart(discCtx, {
            type: 'doughnut',
            data: {
                labels: ['Patent-Free', 'Patented', 'Unknown'],
                datasets: [{
                    data: [104, 10, 5],
                    backgroundColor: [
                        'rgba(0, 230, 118, 0.8)',
                        'rgba(255, 23, 68, 0.8)',
                        'rgba(176, 184, 212, 0.8)'
                    ]
                }]
            },
            options: {
                responsive: true,
                plugins: {
                    legend: { position: 'bottom' }
                }
            }
        });
    }
    
    // Timeline Chart
    const timeCtx = document.getElementById('timelineChart');
    if (timeCtx) {
        new Chart(timeCtx, {
            type: 'line',
            data: {
                labels: ['1985', '1995', '2005', '2015', '2020', '2024'],
                datasets: [{
                    label: 'Research Articles',
                    data: [2, 3, 5, 8, 12, 26],
                    borderColor: 'rgba(0, 168, 232, 1)',
                    backgroundColor: 'rgba(0, 168, 232, 0.1)',
                    fill: true,
                    tension: 0.4
                }]
            },
            options: {
                responsive: true,
                plugins: {
                    legend: { display: false }
                },
                scales: {
                    y: { beginAtZero: true }
                }
            }
        });
    }
}

// Mobile Navigation Toggle
document.getElementById('navToggle')?.addEventListener('click', () => {
    document.querySelector('.nav-menu').classList.toggle('active');
});

