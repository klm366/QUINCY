-- Resolve @sec- and @tbl- references that point at another page.
--
-- Quarto resolves cross-references within one document. Pointing at a label in a
-- different .qmd of a website project renders the literal text "?@sec-x" and warns
-- "Unable to resolve crossref". The numbers do exist, though: Quarto writes them to
-- .quarto/xref while rendering, together with each caption and the file it came from.
--
-- This filter reads that index and rewrites a reference whose target lives on another
-- page into an ordinary link, "Section 3.3" pointing at /that/page.qmd#sec-x, which
-- Quarto then rewrites to the rendered page at the right depth. References within the
-- same page are left untouched, so Quarto keeps numbering those as it always has and
-- there is never a second authority for the same number.
--
-- The index is written during a render, so on a cold build the first render has
-- nothing to read; build-doc.sh renders again, by which point it exists.

local XREF_DIR = nil
local source_file = nil
local index = {}      -- key -> { file = "guide/....qmd", number = "3.3", kind = "sec" }

local prefixes = {
  sec = "Section",
  ch = "Chapter",
  tbl = "Table",
}

-- Whether floats are numbered per chapter ("Table 2.1") or straight through
-- ("Table 3"). Taken from the referring document: the parts of one split document
-- share this setting, and it is the same setting Quarto applies to the references
-- it resolves itself on this page.
local chapters = false

-- [3, 3, 0, 0, ...] -> "3.3"; [2, 0, ...] -> "2"
local function section_number(section)
  local parts = {}
  for _, n in ipairs(section) do
    if n == 0 then break end
    parts[#parts + 1] = tostring(math.floor(n))
  end
  return table.concat(parts, ".")
end

-- Look for a setting in the document, then in the _metadata.yml of each directory
-- up to the project root, then in _quarto.yml: the order Quarto itself merges them.
local function read_setting(path, key)
  local handle = io.open(path, "r")
  if not handle then return nil end
  local content = handle:read("a") or ""
  handle:close()
  local pattern = key:gsub("%-", "%%-")
  local value = content:match("\n%s*" .. pattern .. ":%s*\"([^\"]+)\"")
             or content:match("\n%s*" .. pattern .. ":%s*'([^']+)'")
             or content:match("\n%s*" .. pattern .. ":%s*([^%s\"'#][^\n]*)")
  if value then value = value:gsub("%s+$", "") end
  return value
end

local function find_setting(key)
  if not (quarto and quarto.doc and quarto.doc.input_file and XREF_DIR) then return nil end
  local root = quarto.project.directory
  local value = read_setting(quarto.doc.input_file, key)
  if value then return value end
  local dir = quarto.doc.input_file:match("^(.*)/[^/]*$")
  while dir and #dir >= #root do
    value = read_setting(dir .. "/_metadata.yml", key)
    if value then return value end
    if dir == root then break end
    dir = dir:match("^(.*)/[^/]*$")
  end
  return read_setting(root .. "/_quarto.yml", key)
end

local function read_index()
  if not XREF_DIR then return end
  local index_file = io.open(XREF_DIR .. "/INDEX", "r")
  if not index_file then return end
  local ok, listing = pcall(pandoc.json.decode, index_file:read("a"))
  index_file:close()
  if not ok or type(listing) ~= "table" then return end

  for qmd, outputs in pairs(listing) do
    for _, hash in pairs(outputs) do
      local entry_file = io.open(XREF_DIR .. "/" .. hash, "r")
      if entry_file then
        local parsed, data = pcall(pandoc.json.decode, entry_file:read("a"))
        entry_file:close()
        if parsed and type(data) == "table" and data.entries then
          for _, entry in ipairs(data.entries) do
            local kind = entry.key:match("^(%a+)%-")
            local order = entry.order or {}
            local number
            if kind == "sec" then
              number = section_number(order.section or {})
            elseif chapters and order.section and order.section[1] and order.section[1] > 0 then
              number = math.floor(order.section[1]) .. "." .. math.floor(order.number)
            else
              number = tostring(math.floor(order.number or 0))
            end
            if number and number ~= "" then
              index[entry.key] = { file = qmd, number = number, kind = kind }
            end
          end
        end
      end
    end
  end
end

function Meta(meta)
  if quarto and quarto.project and quarto.project.directory
     and quarto.doc and quarto.doc.input_file then
    local root = quarto.project.directory
    XREF_DIR = root .. "/.quarto/xref"
    source_file = quarto.doc.input_file:gsub("^" .. root:gsub("([^%w])", "%%%1") .. "/?", "")
  end

  if meta.crossref and meta.crossref.chapters then chapters = true end

  -- Follow whatever this document calls its cross-references, so a rewritten one
  -- reads the same as the ones Quarto resolves itself. The prefixes have to be read
  -- from the configuration files: Quarto consumes the `language` block before user
  -- filters run, and what reaches a filter as `labels` is the title block's wording,
  -- not the cross-reference prefixes.
  for short, key in pairs({sec = "crossref-sec-prefix", ch = "crossref-ch-prefix",
                           tbl = "crossref-tbl-prefix"}) do
    local configured = find_setting(key)
    if configured then prefixes[short] = configured end
  end

  read_index()
  return meta
end

function Cite(el)
  -- Only a reference standing on its own; a grouped citation is left alone.
  if #el.citations ~= 1 then return nil end
  local key = el.citations[1].id
  if not key:match("^sec%-") and not key:match("^tbl%-") then return nil end

  local entry = index[key]
  if not entry then return nil end          -- unknown: let Quarto report it
  if entry.file == source_file then return nil end  -- same page: Quarto's job

  local prefix = prefixes[entry.kind] or ""
  if entry.kind == "sec" and not entry.number:match("%.") then
    -- a whole chapter rather than a section within one
    prefix = prefixes.ch
  end

  local text = prefix .. "\u{00A0}" .. entry.number
  return pandoc.Link({pandoc.Str(text)}, "/" .. entry.file .. "#" .. key,
                     "", pandoc.Attr("", {"quarto-xref"}))
end

return {
  { Meta = Meta },
  { Cite = Cite },
}
