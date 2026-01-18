/**
 * PharmaSight™ - Three.js Molecular Animations
 * Advanced molecular visualizations including:
 * - Molecule dissociation/reassembly
 * - Protein-ligand docking
 * - Neurotransmitter synapse animations
 * - Receptor cascade effects
 */

class MolecularAnimationSystem {
    constructor(canvasId) {
        this.canvas = document.getElementById(canvasId);
        if (!this.canvas) return;

        this.scene = null;
        this.camera = null;
        this.renderer = null;
        this.molecules = [];
        this.animationState = 'idle';
        this.currentAnimation = null;

        this.init();
    }

    init() {
        // Setup Three.js scene
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0x0a0a0f);

        // Setup camera
        const aspect = this.canvas.offsetWidth / this.canvas.offsetHeight;
        this.camera = new THREE.PerspectiveCamera(75, aspect, 0.1, 1000);
        this.camera.position.z = 15;

        // Setup renderer
        this.renderer = new THREE.WebGLRenderer({
            canvas: this.canvas,
            antialias: true,
            alpha: true
        });
        this.renderer.setSize(this.canvas.offsetWidth, this.canvas.offsetHeight);
        this.renderer.setPixelRatio(window.devicePixelRatio);

        // Add lights
        this.addLights();

        // Handle resize
        window.addEventListener('resize', () => this.onResize());

        // Start animation loop
        this.animate();
    }

    addLights() {
        // Ambient light
        const ambientLight = new THREE.AmbientLight(0x404040, 2);
        this.scene.add(ambientLight);

        // Directional light
        const directionalLight = new THREE.DirectionalLight(0x00d4ff, 1);
        directionalLight.position.set(5, 5, 5);
        this.scene.add(directionalLight);

        // Point lights for glow effect
        const pointLight1 = new THREE.PointLight(0x00d4ff, 1, 50);
        pointLight1.position.set(-10, 10, 10);
        this.scene.add(pointLight1);

        const pointLight2 = new THREE.PointLight(0x9333ea, 1, 50);
        pointLight2.position.set(10, -10, 10);
        this.scene.add(pointLight2);
    }

    createMolecule(atomPositions, color = 0x00d4ff) {
        const molecule = new THREE.Group();

        // Create atoms
        atomPositions.forEach(pos => {
            const atomGeometry = new THREE.SphereGeometry(0.5, 32, 32);
            const atomMaterial = new THREE.MeshPhongMaterial({
                color: color,
                emissive: color,
                emissiveIntensity: 0.3,
                shininess: 100
            });
            const atom = new THREE.Mesh(atomGeometry, atomMaterial);
            atom.position.set(pos.x, pos.y, pos.z);
            molecule.add(atom);
        });

        // Create bonds
        for (let i = 0; i < atomPositions.length - 1; i++) {
            const start = atomPositions[i];
            const end = atomPositions[i + 1];
            const bond = this.createBond(start, end, color);
            molecule.add(bond);
        }

        return molecule;
    }

    createBond(start, end, color) {
        const direction = new THREE.Vector3().subVectors(
            new THREE.Vector3(end.x, end.y, end.z),
            new THREE.Vector3(start.x, start.y, start.z)
        );
        const length = direction.length();
        const bondGeometry = new THREE.CylinderGeometry(0.1, 0.1, length, 8);
        const bondMaterial = new THREE.MeshPhongMaterial({
            color: color,
            emissive: color,
            emissiveIntensity: 0.2
        });
        const bond = new THREE.Mesh(bondGeometry, bondMaterial);

        bond.position.set(
            (start.x + end.x) / 2,
            (start.y + end.y) / 2,
            (start.z + end.z) / 2
        );

        const axis = new THREE.Vector3(0, 1, 0);
        bond.quaternion.setFromUnitVectors(axis, direction.normalize());

        return bond;
    }

    // ===================================
    // Animation: Molecule Dissociation
    // ===================================
    startDissociationAnimation() {
        this.animationState = 'dissociating';

        // Create a simple molecule
        const atomPositions = [
            { x: 0, y: 0, z: 0 },
            { x: 2, y: 0, z: 0 },
            { x: 4, y: 0, z: 0 },
            { x: 6, y: 0, z: 0 }
        ];

        const molecule = this.createMolecule(atomPositions, 0x00d4ff);
        molecule.position.set(-3, 0, 0);
        this.scene.add(molecule);
        this.molecules.push(molecule);

        // Animate dissociation
        let time = 0;
        const dissociate = () => {
            if (this.animationState !== 'dissociating') return;

            time += 0.02;

            molecule.children.forEach((child, index) => {
                if (child.geometry && child.geometry.type === 'SphereGeometry') {
                    const atomIndex = Math.floor(index / 2);
                    child.position.x += Math.sin(time) * 0.1 * atomIndex;
                    child.position.y += Math.cos(time) * 0.1 * atomIndex;
                }
            });

            molecule.rotation.y += 0.01;

            if (time < 10) {
                requestAnimationFrame(dissociate);
            } else {
                this.startReassemblyAnimation();
            }
        };

        dissociate();
    }

    // ===================================
    // Animation: Molecule Reassembly
    // ===================================
    startReassemblyAnimation() {
        this.animationState = 'reassembling';

        // Create new molecule structure
        const newAtomPositions = [
            { x: 0, y: 0, z: 0 },
            { x: 1.5, y: 1.5, z: 0 },
            { x: 3, y: 0, z: 0 },
            { x: 1.5, y: -1.5, z: 0 }
        ];

        const newMolecule = this.createMolecule(newAtomPositions, 0x9333ea);
        newMolecule.position.set(3, 0, 0);
        newMolecule.scale.set(0, 0, 0);
        this.scene.add(newMolecule);

        // Animate reassembly
        let time = 0;
        const reassemble = () => {
            if (this.animationState !== 'reassembling') return;

            time += 0.05;
            newMolecule.scale.set(time, time, time);
            newMolecule.rotation.y += 0.02;

            if (time < 1) {
                requestAnimationFrame(reassemble);
            } else {
                this.animationState = 'rotating';
                this.startRotationAnimation();
            }
        };

        reassemble();
    }

    // ===================================
    // Animation: Protein-Ligand Docking
    // ===================================
    startDockingAnimation() {
        this.clearScene();
        this.animationState = 'docking';

        // Create protein (larger structure)
        const proteinAtoms = [
            { x: -3, y: 0, z: 0 },
            { x: -1, y: 2, z: 0 },
            { x: 1, y: 2, z: 0 },
            { x: 3, y: 0, z: 0 },
            { x: 1, y: -2, z: 0 },
            { x: -1, y: -2, z: 0 }
        ];
        const protein = this.createMolecule(proteinAtoms, 0x00d4ff);
        protein.position.set(-5, 0, 0);
        protein.scale.set(1.5, 1.5, 1.5);
        this.scene.add(protein);

        // Create ligand (smaller structure)
        const ligandAtoms = [
            { x: 0, y: 0, z: 0 },
            { x: 1, y: 1, z: 0 },
            { x: 2, y: 0, z: 0 }
        ];
        const ligand = this.createMolecule(ligandAtoms, 0xff006e);
        ligand.position.set(8, 0, 0);
        this.scene.add(ligand);

        // Animate docking
        let time = 0;
        const dock = () => {
            if (this.animationState !== 'docking') return;

            time += 0.02;

            // Move ligand towards protein
            ligand.position.x -= 0.1;
            ligand.position.y = Math.sin(time * 2) * 0.5;
            ligand.rotation.z += 0.02;

            // Rotate protein
            protein.rotation.y += 0.005;

            if (ligand.position.x > -2) {
                requestAnimationFrame(dock);
            } else {
                this.animationState = 'docked';
            }
        };

        dock();
    }

    // ===================================
    // Animation: Neurotransmitter Synapse
    // ===================================
    startSynapseAnimation() {
        this.clearScene();
        this.animationState = 'synapse';

        // Create presynaptic neuron (simplified)
        const presynaptic = new THREE.Mesh(
            new THREE.BoxGeometry(4, 2, 1),
            new THREE.MeshPhongMaterial({ color: 0x00d4ff, emissive: 0x00d4ff, emissiveIntensity: 0.2 })
        );
        presynaptic.position.set(-5, 3, 0);
        this.scene.add(presynaptic);

        // Create postsynaptic neuron
        const postsynaptic = new THREE.Mesh(
            new THREE.BoxGeometry(4, 2, 1),
            new THREE.MeshPhongMaterial({ color: 0x00d4ff, emissive: 0x00d4ff, emissiveIntensity: 0.2 })
        );
        postsynaptic.position.set(-5, -3, 0);
        this.scene.add(postsynaptic);

        // Create synaptic cleft
        const cleftGeometry = new THREE.PlaneGeometry(6, 2);
        const cleftMaterial = new THREE.MeshBasicMaterial({
            color: 0x1a1a2e,
            transparent: true,
            opacity: 0.3
        });
        const cleft = new THREE.Mesh(cleftGeometry, cleftMaterial);
        cleft.position.set(-5, 0, -0.5);
        this.scene.add(cleft);

        // Create serotonin molecules
        const serotonins = [];
        for (let i = 0; i < 10; i++) {
            const serotonin = new THREE.Mesh(
                new THREE.SphereGeometry(0.2, 16, 16),
                new THREE.MeshPhongMaterial({ color: 0xff006e, emissive: 0xff006e, emissiveIntensity: 0.5 })
            );
            serotonin.position.set(-5 + Math.random() * 2, 2, Math.random() * 0.5);
            serotonin.userData.velocity = { y: -0.05 - Math.random() * 0.03 };
            this.scene.add(serotonin);
            serotonins.push(serotonin);
        }

        // Create SERT transporters (reuptake)
        const transporters = [];
        for (let i = 0; i < 3; i++) {
            const transporter = new THREE.Mesh(
                new THREE.CylinderGeometry(0.15, 0.15, 0.5, 16),
                new THREE.MeshPhongMaterial({ color: 0x00ff88, emissive: 0x00ff88, emissiveIntensity: 0.3 })
            );
            transporter.position.set(-6 + i * 1, -2, 0);
            transporter.rotation.x = Math.PI / 2;
            this.scene.add(transporter);
            transporters.push(transporter);
        }

        // Animate synapse
        const animate = () => {
            if (this.animationState !== 'synapse') return;

            serotonins.forEach(serotonin => {
                serotonin.position.y += serotonin.userData.velocity.y;

                // Reset if reached bottom
                if (serotonin.position.y < -2) {
                    serotonin.position.y = 2;
                    serotonin.position.x = -5 + Math.random() * 2;
                }

                // Reuptake animation near transporters
                transporters.forEach(transporter => {
                    const distance = serotonin.position.distanceTo(transporter.position);
                    if (distance < 0.5) {
                        serotonin.position.y = 2; // Reset to top
                        serotonin.position.x = -5 + Math.random() * 2;
                    }
                });
            });

            // Pulse transporters
            transporters.forEach((transporter, i) => {
                transporter.scale.y = 1 + Math.sin(Date.now() * 0.005 + i) * 0.1;
            });

            requestAnimationFrame(animate);
        };

        animate();
    }

    // ===================================
    // Animation: Receptor Cascade Effects
    // ===================================
    startReceptorCascadeAnimation() {
        this.clearScene();
        this.animationState = 'cascade';

        // Create receptor
        const receptor = new THREE.Group();

        // Extracellular domain
        const extracellular = new THREE.Mesh(
            new THREE.SphereGeometry(0.8, 16, 16),
            new THREE.MeshPhongMaterial({ color: 0x00d4ff, emissive: 0x00d4ff, emissiveIntensity: 0.3 })
        );
        extracellular.position.y = 2;
        receptor.add(extracellular);

        // Transmembrane domain
        const transmembrane = new THREE.Mesh(
            new THREE.CylinderGeometry(0.3, 0.3, 4, 16),
            new THREE.MeshPhongMaterial({ color: 0x9333ea, emissive: 0x9333ea, emissiveIntensity: 0.2 })
        );
        receptor.add(transmembrane);

        // Intracellular domain
        const intracellular = new THREE.Mesh(
            new THREE.SphereGeometry(0.6, 16, 16),
            new THREE.MeshPhongMaterial({ color: 0x00ff88, emissive: 0x00ff88, emissiveIntensity: 0.3 })
        );
        intracellular.position.y = -2;
        receptor.add(intracellular);

        this.scene.add(receptor);

        // Create ligand
        const ligand = new THREE.Mesh(
            new THREE.SphereGeometry(0.4, 16, 16),
            new THREE.MeshPhongMaterial({ color: 0xff006e, emissive: 0xff006e, emissiveIntensity: 0.5 })
        );
        ligand.position.set(0, 5, 0);
        this.scene.add(ligand);

        // Create G-protein cascade particles
        const gProteins = [];

        // Animate cascade
        let time = 0;
        let cascadeStarted = false;

        const animate = () => {
            if (this.animationState !== 'cascade') return;

            time += 0.02;

            // Ligand binding
            if (ligand.position.y > 2) {
                ligand.position.y -= 0.05;
            } else if (!cascadeStarted) {
                cascadeStarted = true;

                // Trigger cascade - create G-proteins
                for (let i = 0; i < 5; i++) {
                    setTimeout(() => {
                        const gProtein = new THREE.Mesh(
                            new THREE.SphereGeometry(0.2, 16, 16),
                            new THREE.MeshPhongMaterial({ color: 0x00ff88, emissive: 0x00ff88, emissiveIntensity: 0.5 })
                        );
                        gProtein.position.set(0, -2, 0);
                        gProtein.userData.angle = (Math.PI * 2 * i) / 5;
                        gProtein.userData.radius = 0;
                        this.scene.add(gProtein);
                        gProteins.push(gProtein);
                    }, i * 200);
                }
            }

            // Animate G-proteins spreading
            gProteins.forEach(gProtein => {
                gProtein.userData.radius += 0.05;
                const angle = gProtein.userData.angle;
                const radius = gProtein.userData.radius;
                gProtein.position.x = Math.cos(angle) * radius;
                gProtein.position.z = Math.sin(angle) * radius;
                gProtein.position.y = -2 - radius * 0.2;

                // Fade out
                if (gProtein.material.opacity > 0) {
                    gProtein.material.transparent = true;
                    gProtein.material.opacity = Math.max(0, 1 - radius / 5);
                }
            });

            // Rotate receptor
            receptor.rotation.y += 0.01;

            requestAnimationFrame(animate);
        };

        animate();
    }

    // ===================================
    // Utility Methods
    // ===================================
    clearScene() {
        while (this.scene.children.length > 3) { // Keep lights
            const object = this.scene.children[3];
            if (object.geometry) object.geometry.dispose();
            if (object.material) object.material.dispose();
            this.scene.remove(object);
        }
        this.molecules = [];
    }

    startRotationAnimation() {
        const rotate = () => {
            if (this.animationState !== 'rotating') return;

            this.molecules.forEach(molecule => {
                molecule.rotation.x += 0.01;
                molecule.rotation.y += 0.01;
            });

            requestAnimationFrame(rotate);
        };
        rotate();
    }

    onResize() {
        const width = this.canvas.offsetWidth;
        const height = this.canvas.offsetHeight;

        this.camera.aspect = width / height;
        this.camera.updateProjectionMatrix();

        this.renderer.setSize(width, height);
    }

    animate() {
        requestAnimationFrame(() => this.animate());

        // Slowly rotate camera
        this.camera.position.x = Math.sin(Date.now() * 0.0001) * 15;
        this.camera.lookAt(this.scene.position);

        this.renderer.render(this.scene, this.camera);
    }
}

// Initialize when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    // Wait for Three.js to load
    if (typeof THREE !== 'undefined') {
        const molecularSystem = new MolecularAnimationSystem('molecular-canvas');

        // Play button interaction
        const playBtn = document.getElementById('play-animation');
        if (playBtn) {
            let currentAnimationIndex = 0;
            const animations = [
                () => molecularSystem.startDissociationAnimation(),
                () => molecularSystem.startDockingAnimation(),
                () => molecularSystem.startSynapseAnimation(),
                () => molecularSystem.startReceptorCascadeAnimation()
            ];

            playBtn.addEventListener('click', () => {
                animations[currentAnimationIndex]();
                currentAnimationIndex = (currentAnimationIndex + 1) % animations.length;
            });

            // Start with first animation
            molecularSystem.startDissociationAnimation();
        }
    }
});
