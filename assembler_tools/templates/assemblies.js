import {dotplot, assembliesToPaf, loadPaf} from './dotplot.js';

const metadata = fetch('assemblies.json').then(response => response.json()).then(assemblies => {
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

const toggleAllContigGroups = function (assemblyId) {
    // Determine the state of the "select all" checkbox
    const assemblyElement = document.getElementById(`cg-all-${assemblyId}`)
    // Select all contig checkboxes related to the assembly

    const contigs = assemblyElement.parentElement.parentElement.querySelectorAll('.contig');

    // Set the checked state of each contig checkbox based on the "select all" checkbox
    contigs.forEach(contig => {
        const contigName = contig.getAttribute('data-cg')
        toggleContigGroup(contigName, !assemblyElement.checked)
    });
}

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
function aniMatrixGetDendrogramInfo() {
    function extractData(elementId) {
        return JSON.parse(document.getElementById(elementId).textContent)
    }

    function extractLabels(axisSelector) {
        const labels = document.querySelectorAll(axisSelector);
        return Array.from(labels).map(label => label.textContent);
    }

    function isIdentical(textArray1, textArray2) {
        if (textArray1.length !== textArray2.length) return false
        for (let i = 0; i < textArray1.length; i++) {
            if (textArray1[i] !== textArray2[i]) return false
        }
        return true;
    }

    const data = extractData('ani-matrix-data')

    const labels1 = extractLabels('#matplotlib\\.axis_1 text')
    const labels2 = extractLabels('#matplotlib\\.axis_2 text')
    if (!isIdentical(labels1, labels2)) {
        console.error('Dendrogram label mismatch detected. The labels on both axes should be identical for accurate ' +
            'representation. Please verify the input data or the label extraction logic.', {
            axis1Labels: labels1,
            axis2Labels: labels2,
            isIdentical: isIdentical(labels1, labels2)
        })
        alert('The labels on the dendrogram are not identical!')
    }

    return [labels1, data]
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


function createContigGroupContent(contigGroup) {
    const md = window.dataset.contigGroups[contigGroup]  // short for metadata
    return `
    <div class="card">
        <div class="card-header">${md.id}</div>
        <ul class="list-group list-group-flush">
            <li class="list-group-item"><strong>Length</strong>: ${humanBP(md.len)}</li>
            <li class="list-group-item"><strong>GC content</strong>: ${formatAsPercentage(md.gc_rel)}</li>
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

function loadDotplot(event) {
    // Create a div container to hold the generated content
    const div = document.createElement('div');

    // Define the HTML content as a string
    const htmlContent = `
    <div id="dotplot-overlay">
        <button id="dotplot-quit" class="btn btn-danger">Quit</button>
        <div id="dotplot" style="resize: both; overflow: auto; margin: 10px; width: 100%; aspect-ratio: 1 / 1"></div>
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
        console.log({ref, qry})

        const res = assembliesToPaf(ref, qry, 'data')
        res.then(([paf, err, cmd]) => {
            const table = loadPaf(paf)
            dotplot(table, document.getElementById('dotplot'), {title: "Dotplot from minimap2"})
        })
    })
}

/* Show popover on hover */
function aniMatrixInitPopover(dendrogramLabels, dendrogramData) {
    const paths = document.querySelectorAll('#QuadMesh_1 path');
    const popoverMap = new Map();

    function showPopover(popover, rightClick = false) {
        // console.log('showing', {persistent: popover.xx_is_persistent, shown: popover.xx_is_shown})
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
        // console.log('hiding', {persistent: popover.xx_is_persistent, shown: popover.xx_is_shown})
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
            content += createContigGroupContent(labelCol)
        } else {
            content += `Similarity: <strong>${similarity}</strong><br>`
            content += createContigGroupContent(labelCol) + '<br>'
            content += createContigGroupContent(labelRow) + '<br>'
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
            placement: 'auto'
        })

        path.addEventListener('shown.bs.popover', function () {
            document.querySelectorAll('.dotplot-button').forEach(button => {
                button.removeEventListener('click', loadDotplot); // Ensure no duplicate listeners
                button.addEventListener('click', loadDotplot);
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


function initGfavizPopover() {
    document.querySelectorAll('.gfaviz-container').forEach((gfavizContainer) => {
        const assembly = gfavizContainer.getAttribute('data-assembly');
        gfavizContainer.querySelectorAll('svg text').forEach((textElement) => {
            // move them to the front
            textElement.parentNode.parentNode.appendChild(textElement.parentNode);

            // get contig name
            const contigOriginalName = textElement.textContent
            const contigId = `${assembly}@${contigOriginalName}`
            const contigGroupId = window.dataset.contigs[contigId].contig_group

            // set class contig and contig-<contig>
            textElement.classList.add('contig-group');
            textElement.setAttribute('data-cg', contigGroupId)

            // add bootstrap.Popover
            const popoverInstance = new bootstrap.Popover(textElement, {
                trigger: 'manual',
                html: true,
                title: contigGroupId,
                content: 'content',
                container: 'body',
                placement: 'top'
            });

            // Show popover on mouseenter
            textElement.addEventListener('mouseenter', function () {
                const popover = popoverInstance
                if (popover) popover.show()
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


function aniMatrixDeactivateAnnotations() {
    document.getElementById('ani-clustermap-container')
        .querySelectorAll('#axes_3 > [id^="text_"]')
        .forEach((textElement) => {
            textElement.style.pointerEvents = 'none';
        });
}


document.addEventListener('DOMContentLoaded', function () {
    const [dendrogramLabels, dendrogramData] = aniMatrixGetDendrogramInfo();
    aniMatrixDeactivateAnnotations()

    metadata.then(() => {
        aniMatrixInitPopover(dendrogramLabels, dendrogramData);
        initGfavizPopover();
    });
});