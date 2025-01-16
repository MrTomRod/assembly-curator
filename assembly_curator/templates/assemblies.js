import {dotplot, assembliesToPafMinimap, loadPaf} from './dotplot.js';

const sample = document.getElementById('sample').textContent
const metadata = fetch('assembly-curator/assemblies.json').then(response => response.json()).then(assemblies => {
    const contigs = {}
    const contigGroups = {}
    // assemblies is a dictionary of assemblies {'flye': ..., 'canu': ...}
    Object.entries(assemblies).forEach(([assembly, data]) => {
        Object.entries(data.contig_groups).forEach(([contigGroup, data]) => {
            contigGroups[contigGroup] = data
            Object.entries(data.contigs).forEach(([contig, data]) => {
                contigs[contig] = data
            })
        })
    })
    // store as global variable
    window.dataset = {assemblies, contigGroups, contigs}
})

class NoModalError extends Error {
    constructor(message) {
        super(message);
        this.name = 'NoModalError';
    }
}

function showModal(title, bodyHTML, listOfButtons = [["<button type=\"button\" class=\"btn btn-secondary\" data-bs-dismiss=\"modal\">Close</button>", undefined]]) {
    // Create the modal HTML structure
    const modalHTML = `
    <div class="modal fade" tabindex="-1" role="dialog">
        <div class="modal-dialog" role="document">
            <div class="modal-content">
                <div class="modal-header">
                    <h5 class="modal-title">${title}</h5>
                    <button type="button" class="btn-close" data-bs-dismiss="modal" aria-label="Close"></button>
                </div>
                <div class="modal-body">
                    ${bodyHTML}
                </div>
                <div class="modal-footer">
                    ${listOfButtons.map(button => button[0]).join('')}
                </div>
            </div>
        </div>
    </div>`;

    // Append the modal to the body
    const modalElement = document.createElement('div');
    modalElement.innerHTML = modalHTML;
    document.body.appendChild(modalElement);

    // Initialize the modal
    const modal = new bootstrap.Modal(modalElement.querySelector('.modal'));

    // Attach event listeners to buttons
    listOfButtons.forEach(([buttonHTML, clickHandler], index) => {
        if (clickHandler) {
            modalElement.querySelectorAll('.modal-footer button')[index].addEventListener('click', clickHandler);
        }
    });

    // Show the modal
    modal.show();
}

const toggleAllContigGroups = function (assemblyId) {
    // Determine the state of the "select all" checkbox
    const assemblyElement = document.getElementById(`cg-all-${assemblyId}`)

    // Select all contig checkboxes related to the assembly
    const contigGroups = assemblyElement.parentElement.parentElement
        .querySelectorAll('.contig-group');

    // Set the checked state of each contig checkbox based on the "select all" checkbox
    contigGroups.forEach(contigGroup => {
        toggleContigGroup(contigGroup.getAttribute('data-cg'), !assemblyElement.checked)
    });
}
document.querySelectorAll('.assembly-toggler').forEach((checkbox) => {
    checkbox.addEventListener('change', function () {
        toggleAllContigGroups(this.getAttribute('data-assembly'))
    })
})

const toggleContigGroup = function (eventOrElement, isSelected = null) {
    const cg = typeof eventOrElement === 'string' ? eventOrElement : this.getAttribute('data-cg')

    // get all elements that belong to the same contig
    const cgs = Array.from(document.querySelectorAll(`[data-cg="${cg}"]`))

    // if isSelected is not provided, toggle the selection
    if (isSelected === null) {
        isSelected = cgs[0].classList.contains('selected')
    }

    if (isSelected) {
        // select
        cgs.forEach(cg => {
            cg.classList.remove('selected')
            // move selected path elements to the front so that the border is visible
            if (cg.tagName === 'path') cg.parentNode.prepend(cg)
        })
    } else {
        // deselect
        cgs.forEach(cg => {
            cg.classList.add('selected')
            // move deselected path elements to the back so that the border of other paths may be visible
            if (cg.tagName === 'path') cg.parentNode.appendChild(cg)
        })
    }
}

/* Get dendrogram labels */
function aniMatrixGetDendrogramInfo(retryCount = 10) {
    function extractData(elementId) {
        return JSON.parse(document.getElementById(elementId).textContent)
    }

    function extractLabels(labelElements) {
        return Array.from(labelElements).map(label => label.textContent);
    }

    function isIdentical(textArray1, textArray2) {
        if (textArray1.length !== textArray2.length) return false
        for (let i = 0; i < textArray1.length; i++) {
            if (textArray1[i] !== textArray2[i]) return false
        }
        return true;
    }

    const data = extractData('ani-matrix-data')

    const labels1 = extractLabels(document.querySelector('#ani-clustermap-container #matplotlib\\.axis_5').querySelectorAll('text'))
    const labels2 = extractLabels(document.querySelector('#ani-clustermap-container #matplotlib\\.axis_6').querySelectorAll('text'))
    if (!isIdentical(labels1, labels2)) {
        const debugInfo = {axis1Labels: labels1, axis2Labels: labels2, isIdentical: false}
        if (retryCount > 0) {
            console.warn(`Label mismatch detected. Retrying... (${retryCount} retries left)`, debugInfo);
            setTimeout(() => aniMatrixGetDendrogramInfo(retryCount - 1), 500);
            return;
        } else {
            console.error('Dendrogram label mismatch detected after retries.', debugInfo);
            alert('The labels on the dendrogram are not identical!');
        }
    } else {
        return [labels1, data]
    }
}

function formatAsPercentage(floatNumber, decimalPlaces = 2) {
    return `${(floatNumber * 100).toFixed(decimalPlaces)}%`;
}

function humanBP(bp) {
    const magnitude = Math.floor(Math.log10(bp) / 3);
    const shortened = bp / Math.pow(1000, magnitude);
    const units = ['', 'kbp', 'mbp', 'gbp', 'tbp', 'pbp'];
    const unit = units[magnitude] || '';
    return `${shortened.toFixed(1)}${unit}`;
}


function createContigGroupContent(contigGroup, header) {
    const md = window.dataset.contigGroups[contigGroup]  // short for metadata
    header = header || md.id
    const numberOfContigs = Object.keys(md.contigs).length;
    let nOrT, nOrTVal
    if (numberOfContigs === 1) {
        [nOrT, nOrTVal] = ['Topology', md.topology_or_n_contigs]
    } else {
        [nOrT, nOrTVal] = ['Number of contigs', md.contigs.length]
    }
    // loop over md.contigs and get contig.coverage
    const coverage = Object.values(md.contigs).reduce((acc, contig) => acc + contig.coverage, 0)
    return `
    <div class="card">
        <div class="card-header">${header}</div>
        <ul class="list-group list-group-flush">
            <li class="list-group-item"><strong>${nOrT}</strong>: ${nOrTVal}</li>
            <li class="list-group-item"><strong>Length</strong>: ${humanBP(md.len)}</li>
            <li class="list-group-item"><strong>GC content</strong>: ${formatAsPercentage(md.gc_rel)}</li>
            <li class="list-group-item"><strong>Coverage</strong>: ${(coverage)}x</li>
        </ul>
    </div>`
}

function getFasta(contigGroup) {
    const contigGroupRef = window.dataset.contigGroups[contigGroup]
    const assemblerRef = window.dataset.assemblies[contigGroupRef.assembler]
    const pathToFastaRef = assemblerRef.assembly_dir + '/' + assemblerRef.assembly

    let contigs = Object.keys(window.dataset.contigGroups[contigGroup].contigs)

    contigs = contigs.map(contig => contig.substring(contig.indexOf('@') + 1))

    const fastaPromise = fetch(pathToFastaRef).then(response => response.text()).then(fasta => {
        const foundContigs = []
        let resultFasta = ''
        let keep = false
        for (let line of fasta.split('\n')) {
            if (line.startsWith('>')) {
                line = line.split(/\s+/)[0]  // split by whitespace and take the first part
                if (contigs.includes(line.slice(1))) {
                    foundContigs.push(line.slice(1))
                    keep = true
                    resultFasta += line + '\n'
                } else {
                    keep = false
                }
            } else {
                if (keep) {
                    resultFasta += line + '\n'
                }
            }
        }
        // check if all contigs were found
        if (foundContigs.length !== contigs.length) {
            const msg = `Not all contigs were found in the fasta file. Missing: ${contigs.filter(contig => !foundContigs.includes(contig))}`
            alert(msg)
            throw new Error(msg)
        }
        return resultFasta
    })
    return fastaPromise
}

function loadDotplot() {
    // Create a div container to hold the generated content
    const div = document.createElement('div');

    // Define the HTML content as a string
    const htmlContent = `
    <div id="dotplot-overlay">
        <div id="dotplot" style=""></div>
        <button id="dotplot-quit" class="btn btn-danger">Quit</button>
    </div>`;

    // Add the div to the document body
    div.innerHTML = htmlContent;
    document.body.appendChild(div);

    // Make the quit button functional
    document.getElementById('dotplot-quit').addEventListener('click', () => {
        const overlay = document.getElementById('dotplot-overlay');
        if (overlay) {
            overlay.remove();
        }
    });

    const fastaRef = getFasta(this.dataset.ref)
    const fastaQry = getFasta(this.dataset.qry)

    Promise.all([fastaRef, fastaQry]).then(([ref, qry]) => {
        const res = assembliesToPafMinimap(ref, qry, 'data')
        res.then(([paf, err, cmd]) => {
            const table = loadPaf(paf)
            dotplot(table, document.getElementById('dotplot'), {title: "Dotplot from minimap2"})
        })
    })
}

/**
 * Extracts a random subsequence of a specified length from a given sequence.
 *
 * @param {string} sequence - The input DNA/RNA sequence.
 * @param {number} [length=1000] - The desired length of the random subsequence (default: 1000 bp).
 * @param {string} newHeader - An optional new FASTA header
 * @returns {string} - A random subsequence of the specified length, or the full sequence if it is shorter than the desired length.
 */
function extractRandomSubsequence(fasta, length = 1000, newHeader = undefined) {
    // Split by ">" to separate contigs, and filter out empty entries
    const contigs = fasta.split(">").filter(entry => entry.trim() !== "")

    // Select a random contig
    const randomContig = contigs[Math.floor(Math.random() * contigs.length)]

    // Split header and sequence
    const [header, ...sequenceLines] = randomContig.split("\n")

    // if newHeader is not provided, use the original header
    if (!newHeader) newHeader = header

    // Remove newlines and whitespace
    let sequence = sequenceLines.join("").replace(/\s+/g, "")

    // If the sequence is shorter than the desired length, return the whole sequence
    if (sequence.length > length) {
        const start = Math.floor(Math.random() * (sequence.length - length + 1))
        sequence = sequence.slice(start, start + length)
    }
    return `>${newHeader}\n${sequence}`
}

function blastFasta() {
    const contigGroup = this.dataset.contigGroup
    const fastaPromise = getFasta(contigGroup)
    const nBases = 1000
    const header = `${sample}: ${contigGroup} (first ${nBases} bp)`
    fastaPromise.then(fasta => {
        fasta = extractRandomSubsequence(fasta, nBases, header)

        /*
        // Simply open the BLAST page with the encoded FASTA sequence:
        const fastaEncoded = encodeURIComponent(fasta)
        window.open(`https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome&QUERY=${fastaEncoded}`)
        */

        // Send the FASTA sequence to the BLAST API:
        fetch("https://blast.ncbi.nlm.nih.gov/Blast.cgi", {
            // https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html
            method: "POST",
            headers: {
                "Content-Type": "application/x-www-form-urlencoded",
            },
            body: new URLSearchParams({
                CMD: "Put",
                QUERY: fasta,
                DATABASE: "core_nt",
                PROGRAM: "blastn",
                MEGABLAST: "on",
            }).toString(),
        }).then(response => response.text()).then(html => {
            // Parse the HTML response
            const parser = new DOMParser();
            const doc = parser.parseFromString(html, "text/html");
            // Select the RID input element
            const ridInput = doc.querySelector('input[name="RID"]');
            if (ridInput) {
                const rid = ridInput.value;
                console.log("Job submitted. RID:", rid);
                // Open the results page
                window.open(`https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=${rid}`, "_blank");
            } else {
                alert("Failed to submit job to NCBI BLAST. See error console for more details.");
                console.error("Failed to submit job to NCBI BLAST.", {ridInput, rid, html});
            }
        })
            .catch(error => console.error("Error:", error));

        /*
        // Send the FASTA sequence to the EBI BLAST API:
        // Problem: Much slower than NCBI
        fetch("https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run/", {
            // https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/parameters
            method: "POST",
            headers: {
                "Content-Type": "application/x-www-form-urlencoded",
            },
            body: new URLSearchParams({
                "sequence": fasta,
                "stype": "dna",
                "database": "em_pro",
                "program": "blastn",
                "task": "megablast",
                "email": "suranj@ebi.ac.uk",
                "title": header
            }).toString()
        })
            .then(response => response.text()).then((jobId) => {
            console.log("Job submitted. Job ID:", jobId)
            window.open(`https://www.ebi.ac.uk/jdispatcher/sss/ncbiblast/summary?jobId=${jobId}`, "_blank");
        })
            .catch(error => console.error("Error:", error))
        */
    })
}

/* Show popover on hover */
function aniMatrixInitPopover(dendrogramLabels, dendrogramData) {
    const paths = document.querySelectorAll('#QuadMesh_3 path');
    const popoverMap = new Map();

    function showPopover(popover, rightClick = false) {
        if (rightClick) {
            popover.xx_is_persistent = !popover.xx_is_persistent
            if (!popover.xx_is_persistent) {
                hidePopover(popover)
                return
            }
        }
        if (popover.xx_is_shown) return
        popover.xx_is_shown = true
        popover.show()
    }

    function hidePopover(popover) {
        if (popover.xx_is_persistent) return
        if (popover.xx_is_shown) {
            popover.xx_is_shown = false
            popover.hide()
        }
    }

    function forceHidePopover(popover) {
        if (popover.xx_is_shown) popover.hide()
        popover.xx_is_persistent = false
        popover.xx_is_shown = false
    }

    const getLabels = function (id) {
        const index = typeof id === 'string' ? parseInt(id.split('-')[1]) : id;
        const labelCol = dendrogramLabels[index % dendrogramLabels.length];
        const labelRow = dendrogramLabels[Math.floor(index / dendrogramLabels.length)];
        return [labelCol, labelRow];
    }

    const titleFunction = function (element) {
        const labelCol = element.getAttribute('data-label-col')
        const labelRow = element.getAttribute('data-label-row')
        if (labelCol === labelRow) {
            return `${labelCol}`
        } else {
            return `${labelCol} x ${labelRow}`
        }
    }

    const contentFunction = function (element) {
        const labelCol = element.getAttribute('data-label-col')
        const labelRow = element.getAttribute('data-label-row')
        const similarity = dendrogramData[labelCol][labelRow]
        let content = ''
        if (labelCol === labelRow) {
            content += createContigGroupContent(labelCol) + '<br>'
            content += `
            <button type="button" class="btn btn-primary blast-button" data-contig-group="${labelCol}">Blast 1000bp</button>`
        } else {
            content += `Similarity: <strong>${similarity}</strong><br>`
            content += createContigGroupContent(labelCol, `Col: ${labelCol}`) + '<br>'
            content += createContigGroupContent(labelRow, `Row: ${labelRow}`) + '<br>'
        }
        // add button to create dotplot
        content += `
            <button type="button" class="btn btn-primary dotplot-button" data-ref="${labelCol}" data-qry="${labelRow}">Dotplot</button>`
        return content
    }

    paths.forEach((path, index) => {
        // Set a unique ID for each path to target with the popover
        const [labelCol, labelRow] = getLabels(index)
        path.setAttribute('id', `path-${index}`);
        path.setAttribute('data-label-col', labelCol);
        path.setAttribute('data-label-row', labelRow);

        // Add click event listener to toggle border
        if (labelCol === labelRow) {
            path.classList.add('contig-group');
            path.setAttribute('data-cg', labelCol);
            path.addEventListener('click', toggleContigGroup);
        }

        const popoverInstance = new bootstrap.Popover(path, {
            trigger: 'manual',
            animation: false,
            html: true,
            sanitizeFn: (content) => content,
            title: titleFunction,
            content: contentFunction,
            container: 'body',
            placement: 'bottom'
        })

        path.addEventListener('shown.bs.popover', function () {
            document.querySelectorAll('.dotplot-button').forEach(button => {
                button.removeEventListener('click', loadDotplot); // Ensure no duplicate listeners
                button.addEventListener('click', loadDotplot);
            });
            document.querySelectorAll('.blast-button').forEach(button => {
                button.removeEventListener('click', blastFasta); // Ensure no duplicate listeners
                button.addEventListener('click', blastFasta);
            });
        });

        // custom variables
        popoverInstance.xx_is_shown = false
        popoverInstance.xx_is_persistent = false

        // add to popoverMap
        popoverMap.set(path, popoverInstance);

        // Show popover on mouseenter
        path.addEventListener('mouseenter', function (event) {
            showPopover(popoverInstance);
        });

        // Hide popover on mouseleave
        path.addEventListener('mouseleave', function (event) {
            hidePopover(popoverInstance);
        });

        // Show persistent popover on right-click
        path.addEventListener('contextmenu', function (event) {
            event.preventDefault();
            showPopover(popoverInstance, true);
        });
    })

    // Hide all persistent popovers if clicked outside
    document.addEventListener('click', function (event) {
        // Ensure it's a left-click
        if (event.button !== 0) return
        // Ensure the click is not on a popover or the dotplot is being used
        const closestPopoverElement = event.target.closest('.popover, #dotplot-overlay')
        if (closestPopoverElement) return
        popoverMap.forEach(forceHidePopover)
    })

    // Hide all persistent popovers if escape key is pressed
    document.addEventListener('keydown', function (event) {
        if (event.key === 'Escape' || event.key === 'Esc') {
            const dotplotOverlay = document.getElementById('dotplot-overlay')
            if (dotplotOverlay) {
                dotplotOverlay.remove()
            } else {
                popoverMap.forEach(forceHidePopover)
            }
        }
    })
}

function makeDraggable(svg) {
    svg.addEventListener('mousedown', startDrag);
    svg.addEventListener('mousemove', drag);
    svg.addEventListener('mouseup', endDrag);
    svg.addEventListener('mouseleave', endDrag);

    let selectedElement, offset;


    function getMousePosition(evt) {
        const CTM = selectedElement.parentElement.getCTM();
        return {
            x: (evt.clientX - CTM.e) / CTM.a, y: (evt.clientY - CTM.f) / CTM.d
        };
    }


    function startDrag(evt) {
        if (evt.target.classList.contains('draggable')) {
            selectedElement = evt.target;
            offset = getMousePosition(evt);
            offset.x -= parseFloat(selectedElement.getAttributeNS(null, "x"));
            offset.y -= parseFloat(selectedElement.getAttributeNS(null, "y"));
        }
    }

    function drag(evt) {
        if (selectedElement) {
            const coord = getMousePosition(evt);
            evt.preventDefault();
            selectedElement.setAttributeNS(null, "x", coord.x - offset.x);
            selectedElement.setAttributeNS(null, "y", coord.y - offset.y);
        }
    }

    function endDrag(evt) {
        selectedElement = null;
    }

    return function isDragging() {
        return selectedElement !== null;
    }
}

function gfavizInitPopover() {
    document.querySelectorAll('.gfaviz-container').forEach((gfavizContainer) => {
        const svg = gfavizContainer.querySelector('svg')
        if (!svg) return

        // add 20 more pixels at the bottom: unset height and width, edit viewport
        svg.removeAttribute('height')
        svg.removeAttribute('width')
        svg.setAttribute('viewBox', `0 0 ${svg.viewBox.baseVal.width} ${svg.viewBox.baseVal.height + 40}`)

        const assembly = gfavizContainer.getAttribute('data-assembly');
        const isDragging = makeDraggable(svg)

        const uniqueTextContents = new Set();

        gfavizContainer.querySelectorAll('svg text').forEach((textElement) => {
            // get contig name
            const contigOriginalName = textElement.textContent
            if (uniqueTextContents.has(contigOriginalName)) {
                textElement.remove(); // The text elements are duplicated. No idea why this happens.
                return
            } else {
                uniqueTextContents.add(contigOriginalName);
            }

            const contigId = `${assembly}@${contigOriginalName}`

            // move them to the front
            textElement.parentNode.parentNode.appendChild(textElement.parentNode);
            textElement.classList.add('draggable');

            // sometimes, a contig in the svg is not in the metadata because it was filtered out during polishing (Flye)
            if (!window.dataset.contigs[contigId]) {
                // add Popover with warning
                const popoverInstance = new bootstrap.Popover(textElement, {
                    trigger: 'hover',
                    title: 'Warning',
                    content: `Contig ${contigOriginalName} is missing in the assembly FASTA file.`,
                    container: 'body',
                    placement: 'top',
                    customClass: 'popover-danger'
                });
                return
            }

            const contigGroupData = window.dataset.contigs[contigId]
            const contigGroupId = contigGroupData.contig_group

            // set stroke to cluster color
            textElement.style.stroke = contigGroupData.cluster_color_rgb

            // set class contig and contig-<contig>
            textElement.classList.add('contig-group');
            textElement.setAttribute('data-cg', contigGroupId)

            // add bootstrap.Popover
            const popoverInstance = new bootstrap.Popover(textElement, {
                trigger: 'manual',
                html: true,
                title: contigGroupId,
                content: createContigGroupContent(contigGroupId),
                container: 'body',
                placement: 'top'
            });

            // Show popover on mouseenter
            textElement.addEventListener('mouseenter', function () {
                const popover = popoverInstance
                if (popover && !isDragging()) popover.show()
            });

            // Hide popover on mouseleave
            textElement.addEventListener('mouseleave', function () {
                const popover = popoverInstance
                if (popover) popover.hide()
            })

            textElement.addEventListener('click', toggleContigGroup)
        });
    });
};

function dotplotsInitPopover() {
    document.querySelectorAll('#cluster-tabs-content svg').forEach((svg) => {
        svg.querySelectorAll('[id^="dotplot - "]').forEach((path) => {
            const [_, labelCol, labelRow] = path.id.split(' - ')
            path.setAttribute('data-label-col', labelCol);
            path.setAttribute('data-label-row', labelRow);

            if (labelCol === labelRow) {
                path.classList.add('contig-group');
                path.setAttribute('data-cg', labelCol);
                path.addEventListener('click', toggleContigGroup);
            }

        })
    })
}


// Trigger toggleContigGroup on click on contigs in the table
function toggleContigGroupTable() {
    document.querySelectorAll('#row-contigs [data-cg]').forEach((contig) => {
        contig.addEventListener('click', toggleContigGroup);
    })
}


function aniMatrixDeactivateAnnotations() {
    document.getElementById('ani-clustermap-container')
        .querySelectorAll('#axes_5 > [id^="text_"]')
        .forEach((textElement) => {
            textElement.style.pointerEvents = 'none';
        });
}

// click on #btn-curate will send the selected contigs to the server (/curate)
document.getElementById('btn-export').addEventListener('click', function () {
    const payload = {
        sample: sample,
        contigs: Array.from(document.querySelectorAll('#row-contigs .contig-group.selected')).map(contig => contig.getAttribute('data-cg')),
    }

    if (payload.contigs.length === 0) {
        alert('No contigs selected!')
        return
    }

    // Add iframe to #curate-div
    // Similar to <iframe src="/curate?sample=${payload.sample}&contigs=${payload.contigs.join(',')}" title="Curate"></iframe>
    const targetDiv = document.getElementById('export-div')
    const iframe = document.createElement('iframe')

    const encodedSample = encodeURIComponent(payload.sample);
    const encodedContigs = payload.contigs.map(contig => encodeURIComponent(contig)).join(',');
    iframe.src = `/export?sample=${encodedSample}&contig_groups=${encodedContigs}`;
    iframe.title = `Export ${payload.sample}`
    targetDiv.innerHTML = ''
    targetDiv.appendChild(iframe)
})

function fetchAndReplace(imgElement) {
    const src = imgElement.getAttribute('src');

    if (!imgElement.complete || imgElement.naturalWidth === 0) {
        console.info(`Not replacing ${src} with SVG because it does not exist.`);
        return
    }
    /* Replace <img src=...> with fetched SVG so it can be made interactive */
    return fetch(src).then(response => response.text()).then(svg => {
        const div = document.createElement('div');
        div.innerHTML = svg;
        const svgElement = div.querySelector('svg');
        imgElement.parentNode.replaceChild(svgElement, imgElement);
    }).catch(error => {
        console.error(`Error fetching SVG ${src}:`, error);
    })
}


function toggleDotPlotTabs() {
    /* Hide tab content if clicked a second time */
    let previouslySelected = undefined
    document.querySelectorAll('#cluster-tabs .nav-link').forEach(tab => {
        const target = document.querySelector(tab.getAttribute('data-bs-target'));
        /* Hide tab content if clicked a second time */
        tab.addEventListener('click', function () {
            if (previouslySelected === target) {
                target.classList.toggle('active');
            } else {
                previouslySelected = target
            }
        });
    })
}

function calculateATGCgradient(div) {
    // Parse the data-atgc attribute
    const atgcData = JSON.parse(div.getAttribute('data-atgc').replace(/'/g, '"'));
    const [a, t, g, c] = [atgcData.A, atgcData.T, atgcData.G, atgcData.C];
    const total = a + t + g + c;

    // Calculate the percentages
    const [aRel, tRel, gRel, cRel] = [a / total * 100, t / total * 100, g / total * 100, c / total * 100];

    // Calculate cumulative percentages for the conic-gradient
    const [aEnd, tEnd, gEnd, cEnd] = [aRel, aRel + tRel, aRel + tRel + gRel, aRel + tRel + gRel + cRel];


    // Create the linear gradient CSS value
    const gradient = `linear-gradient(
                to right,
                blue 0%,         blue ${aEnd}%,
                yellow ${aEnd}%, yellow ${tEnd}%,
                green ${tEnd}%,  green ${gEnd}%,
                red ${gEnd}%,    red 100%
            )`;
    const tooltipText = `A: ${aRel.toFixed(2)}%<br>T: ${tRel.toFixed(2)}%<br>G: ${gRel.toFixed(2)}%<br>C: ${cRel.toFixed(2)}%`;
    return [gradient, tooltipText];
}

function addATGCgradient(div, tooltip = true) {
    const [gradient, tooltipText] = calculateATGCgradient(div);
    // Apply the gradient as the background of the div
    div.style.setProperty('--line-gradient', gradient);

    if (tooltip) {
        // Set up the tooltip using Bootstrap's data-bs-toggle attribute
        div.setAttribute('data-bs-toggle', 'tooltip');
        div.setAttribute('data-bs-html', 'true');
        div.setAttribute('title', tooltipText);

        // Initialize the tooltip
        new bootstrap.Tooltip(div);
    }
}

function activateAniZoom() {
    const container = document.getElementById('ani-clustermap-container');
    const zoomInButton = document.getElementById('ani-zoom-in');
    const zoomOutButton = document.getElementById('ani-zoom-out');
    let scale = 1; // Scale factor for zooming

    // Function to adjust the SVG width based on scale
    function adjustSVGWidth() {
        const svg = container.querySelector('svg') || container.querySelector('img')
        svg.style.width = `${scale * 100}%`;
    }

    // Zoom In Function
    zoomInButton.addEventListener('click', () => {
        scale += 0.1;
        adjustSVGWidth();
    });

    // Zoom Out Function
    zoomOutButton.addEventListener('click', () => {
        if (scale > 0.1) { // Prevent zooming out too much
            scale -= 0.1;
            adjustSVGWidth();
        }
    });

    // Keyboard event listener for '+' and '-' keys
    document.addEventListener('keydown', (event) => {
        if (event.key === '+' || event.key === '=') { // '+' key
            scale += 0.1;
            adjustSVGWidth();
        } else if (event.key === '-') { // '-' key
            if (scale > 0.1) {
                scale -= 0.1;
                adjustSVGWidth();
            }
        }
    });
}

function extractHeaders(fastaContent) {
    // Split the content by newlines
    const lines = fastaContent.split('\n');

    // Initialize the dictionary to store the headers and their metadata
    const headersDict = {};

    // Filter lines that start with ">"
    const headers = lines.filter(line => line.startsWith('>'));

    // Loop over the headers and extract metadata
    headers.forEach(header => {
        // Extract the identifier (e.g., "FAM2878-p1-1_scf1")
        const identifier = header.split(' ')[0].substring(1);

        // Extract the metadata (e.g., "length=2498046", "topology=circular", ...)
        const metadata = header.match(/\[([^\]]+)\]/g).reduce((acc, item) => {
            const [key, value] = item.slice(1, -1).split('=');
            acc[key] = value;
            return acc;
        }, {});

        // Store the metadata in the dictionary
        headersDict[identifier] = metadata;
    });

    return headersDict;
}

function cssEscape(str) {
    return str.replace(/([ #;?%&,.+*~\':"!^$[\]()=>|\/@])/g, '\\$1');
}

function selectBasedOnHybridFasta() {
    fetch('./assembly-curator/hybrid.fasta')
        .then(response => {
            if (response.ok) {
                console.info('attempting to select contigs based on hybrid.fasta...');
                return response.text();
            } else {
                console.info('hybrid.fasta does not exist.');
                return Promise.reject('File not found');
            }
        })
        .then(fastaContent => {
            const headersDict = extractHeaders(fastaContent)
            const rowContigsElement = document.getElementById('row-contigs');
            Object.entries(headersDict).forEach(([identifier, dataDict]) => {
                const [assembler, contigId] = [dataDict['assembler'], dataDict['old-id']]
                const contigAssId = `${assembler}@${contigId}`
                const matchingElements = rowContigsElement.querySelector(`.${cssEscape(contigAssId)}`)
                if (!matchingElements) {
                    const message = `Contig ${contigAssId} not found in the table.`
                    document.querySelector('.error-container')
                        .innerHTML += `<div class="alert alert-danger" role="alert">${message}</div>`
                    return
                }
                const contigGroup = matchingElements.closest('[data-cg]').getAttribute('data-cg')
                toggleContigGroup(contigGroup, false)
            })
        })
        .catch(error => {
            if (error !== 'File not found') {
                console.error('Error fetching hybrid.fasta:', error);
            }
        });
}

function writeFile(targetFile, content) {
    return fetch(targetFile, {
        method: 'PUT',
        headers: {
            'Content-Type': 'text/plain',
        },
        body: content,
    })
}

function loadNoteMd() {
    let noteContainer = document.querySelector('.note-container');
    if (!noteContainer) {
        const firstH1 = document.querySelector('h1');
        if (firstH1) {
            noteContainer = document.createElement('div');
            noteContainer.classList.add('note-container', 'container');
            firstH1.parentNode.insertBefore(noteContainer, firstH1.nextSibling);
        } else {
            console.error("No <h1> element found to insert the note container.");
            return;
        }
    }
    noteContainer.innerHTML = '';  // Empty the noteContainer

    // Style the note-container to center its content
    noteContainer.style.display = 'flex';
    noteContainer.style.justifyContent = 'center';

    fetch('./note.md')
        .then(function (response) {
            if (!response.ok) {
                const button = document.createElement('button');
                button.classList.add('btn', 'btn-secondary');
                button.title = 'add note.md';
                button.innerHTML = '<i class="bi bi-pencil-square"></i>';
                button.addEventListener('click', function () {
                    writeFile('./note.md', '### Empty note\n').then(() => {
                        loadNoteMd();
                    })
                        .catch(() => {
                            showModal('Error', `<p>Failed to write to <code>./note.md</code>.</p>`);
                        })
                    button.remove();
                });
                document.querySelector('h1>.btn-group').appendChild(button);
                throw new NoModalError('note.md not found');
            }
            return response.text();
        })
        .then(function (markdownText) {
            // Convert the Markdown to HTML using Snarkdown
            const htmlContent = snarkdown(markdownText);

            // Create a Bootstrap card
            const card = document.createElement('div');
            card.classList.add('card', 'text-bg-warning');
            card.style.textAlign = 'left'; // Align text to the left
            card.style.width = 'max-content'; // Set card width to max-content
            card.style.minWidth = '15rem';

            const cardHeader = document.createElement('div');
            cardHeader.classList.add('card-header');
            cardHeader.textContent = 'note.md';

            const cardBody = document.createElement('div');
            cardBody.classList.add('card-body');

            const cardContent = document.createElement('div');
            cardContent.innerHTML = htmlContent;

            const editButton = document.createElement('button');
            editButton.classList.add('btn', 'btn-primary', 'mt-3');
            editButton.textContent = 'Edit';
            editButton.type = 'button';

            // Append cardContent before adding functionality
            cardBody.appendChild(cardContent);

            function createTextarea() {
                const textarea = document.createElement('textarea');
                textarea.classList.add('form-control', 'mt-3'); // Bootstrap styling
                textarea.value = markdownText; // Populate with original Markdown
                textarea.style.resize = 'both';
                textarea.style.minHeight = '10rem';
                cardBody.replaceChild(textarea, cardContent);

                // Change the button text to "Save"
                editButton.textContent = 'Save';

                // Change the button's behavior to save the edited content
                editButton.removeEventListener('click', createTextarea);
                editButton.addEventListener('click', function () {
                    writeFile('./note.md', textarea.value)
                        .then(() => {
                            // Empty noteContainer
                            noteContainer.innerHTML = '';
                            // Reload the note.md file
                            loadNoteMd();
                        })
                        .catch(() => {
                            showModal('Error', `<p>Failed to write to <code>./note.md</code>.</p>`);
                        });
                });
            }

            editButton.addEventListener('click', createTextarea);

            // Assemble the card
            cardBody.appendChild(editButton);
            card.appendChild(cardHeader);
            card.appendChild(cardBody);

            // Append the card to the noteContainer
            noteContainer.appendChild(card);
        })
        .catch(function (error) {
            if (error instanceof NoModalError) {
                console.info(error.message);
            } else {
                showModal('Error', `<p>Failed to load or render <code>./note.md</code>.</p>`);
            }
        });
}

document.addEventListener('DOMContentLoaded', function () {
    try {
        activateAniZoom();
    } catch (error) {
        console.error('activateAniZoom failed, but continuing execution.', error);
    }
    toggleContigGroupTable()
    toggleDotPlotTabs()

    loadNoteMd()

    metadata.then(() => {
        fetchAndReplace(document.getElementById('ani-matrix-svg')).then(() => {
            const [dendrogramLabels, dendrogramData] = aniMatrixGetDendrogramInfo();
            aniMatrixDeactivateAnnotations()
            aniMatrixInitPopover(dendrogramLabels, dendrogramData);
        })

        const replaceGfaviz = Promise.all(Array.from(document.querySelectorAll('.gfaviz-svg')).map(fetchAndReplace)).then(() => {
            gfavizInitPopover();
        });

        const replaceDotplots = Promise.all(Array.from(document.querySelectorAll('.dotplot-svg')).map(fetchAndReplace)).then(() => {
            dotplotsInitPopover()
        })

        Promise.all([replaceGfaviz, replaceDotplots]).then(() => {
            // if hybrid.fasta exists, select the chosen ContigGroups using toggleContigGroup
            selectBasedOnHybridFasta()
        });
    });

    document.querySelectorAll('.assembly-atgc').forEach(addATGCgradient);
    document.querySelectorAll('.atgc-indicator').forEach(addATGCgradient);
});