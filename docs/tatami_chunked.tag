<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.9.8">
  <compound kind="file">
    <name>CustomDenseChunkedMatrix.hpp</name>
    <path>tatami_chunked/</path>
    <filename>CustomDenseChunkedMatrix_8hpp.html</filename>
    <class kind="struct">tatami_chunked::CustomDenseChunkedOptions</class>
    <class kind="class">tatami_chunked::CustomDenseChunkedMatrix</class>
    <namespace>tatami_chunked</namespace>
  </compound>
  <compound kind="file">
    <name>CustomSparseChunkedMatrix.hpp</name>
    <path>tatami_chunked/</path>
    <filename>CustomSparseChunkedMatrix_8hpp.html</filename>
    <class kind="struct">tatami_chunked::CustomSparseChunkedOptions</class>
    <class kind="class">tatami_chunked::CustomSparseChunkedMatrix</class>
    <namespace>tatami_chunked</namespace>
  </compound>
  <compound kind="file">
    <name>LruSlabCache.hpp</name>
    <path>tatami_chunked/</path>
    <filename>LruSlabCache_8hpp.html</filename>
    <class kind="class">tatami_chunked::LruSlabCache</class>
    <namespace>tatami_chunked</namespace>
  </compound>
  <compound kind="file">
    <name>mock_dense_chunk.hpp</name>
    <path>tatami_chunked/</path>
    <filename>mock__dense__chunk_8hpp.html</filename>
    <class kind="struct">tatami_chunked::MockSimpleDenseChunk</class>
    <class kind="struct">tatami_chunked::SimpleDenseChunkWrapper</class>
    <class kind="struct">tatami_chunked::MockSubsetDenseChunk</class>
    <namespace>tatami_chunked</namespace>
  </compound>
  <compound kind="file">
    <name>mock_sparse_chunk.hpp</name>
    <path>tatami_chunked/</path>
    <filename>mock__sparse__chunk_8hpp.html</filename>
    <class kind="struct">tatami_chunked::MockSimpleSparseChunk</class>
    <class kind="struct">tatami_chunked::MockSimpleSparseChunk::Workspace</class>
    <class kind="struct">tatami_chunked::SimpleSparseChunkWrapper</class>
    <class kind="struct">tatami_chunked::MockSubsetSparseChunk</class>
    <class kind="struct">tatami_chunked::MockSubsetSparseChunk::Workspace</class>
    <namespace>tatami_chunked</namespace>
  </compound>
  <compound kind="file">
    <name>OracleSlabCache.hpp</name>
    <path>tatami_chunked/</path>
    <filename>OracleSlabCache_8hpp.html</filename>
    <class kind="class">tatami_chunked::OracleSlabCache</class>
    <namespace>tatami_chunked</namespace>
  </compound>
  <compound kind="file">
    <name>SubsettedOracleSlabCache.hpp</name>
    <path>tatami_chunked/</path>
    <filename>SubsettedOracleSlabCache_8hpp.html</filename>
    <class kind="struct">tatami_chunked::SubsettedOracleSlabSubset</class>
    <class kind="class">tatami_chunked::SubsettedOracleSlabCache</class>
    <namespace>tatami_chunked</namespace>
  </compound>
  <compound kind="file">
    <name>tatami_chunked.hpp</name>
    <path>tatami_chunked/</path>
    <filename>tatami__chunked_8hpp.html</filename>
    <namespace>tatami_chunked</namespace>
  </compound>
  <compound kind="file">
    <name>typical_slab_cache.hpp</name>
    <path>tatami_chunked/</path>
    <filename>typical__slab__cache_8hpp.html</filename>
    <class kind="struct">tatami_chunked::TypicalSlabCacheWorkspace</class>
    <namespace>tatami_chunked</namespace>
  </compound>
  <compound kind="class">
    <name>tatami_chunked::CustomDenseChunkedMatrix</name>
    <filename>classtatami__chunked_1_1CustomDenseChunkedMatrix.html</filename>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>typename Chunk_</templarg>
    <base>Matrix&lt; Value_, Index_ &gt;</base>
    <member kind="function">
      <type></type>
      <name>CustomDenseChunkedMatrix</name>
      <anchorfile>classtatami__chunked_1_1CustomDenseChunkedMatrix.html</anchorfile>
      <anchor>ac60c324ee232ab578353f879e1540083</anchor>
      <arglist>(Index_ mat_nrow, Index_ mat_ncol, Index_ chunk_nrow, Index_ chunk_ncol, std::vector&lt; Chunk_ &gt; chunks, bool row_major, const CustomDenseChunkedOptions &amp;opt)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami_chunked::CustomDenseChunkedOptions</name>
    <filename>structtatami__chunked_1_1CustomDenseChunkedOptions.html</filename>
    <member kind="variable">
      <type>size_t</type>
      <name>maximum_cache_size</name>
      <anchorfile>structtatami__chunked_1_1CustomDenseChunkedOptions.html</anchorfile>
      <anchor>ac867e9ee92d981a4b64bdcd3922ba0ef</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>require_minimum_cache</name>
      <anchorfile>structtatami__chunked_1_1CustomDenseChunkedOptions.html</anchorfile>
      <anchor>a371c078f772cc4b494988a7b96d6798e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami_chunked::CustomSparseChunkedMatrix</name>
    <filename>classtatami__chunked_1_1CustomSparseChunkedMatrix.html</filename>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>typename Chunk_</templarg>
    <base>Matrix&lt; Value_, Index_ &gt;</base>
    <member kind="function">
      <type></type>
      <name>CustomSparseChunkedMatrix</name>
      <anchorfile>classtatami__chunked_1_1CustomSparseChunkedMatrix.html</anchorfile>
      <anchor>aaf8de32d64a64a80f2c8a7d3467871ee</anchor>
      <arglist>(Index_ mat_nrow, Index_ mat_ncol, Index_ chunk_nrow, Index_ chunk_ncol, std::vector&lt; Chunk_ &gt; chunks, bool row_major, const CustomSparseChunkedOptions &amp;opt)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami_chunked::CustomSparseChunkedOptions</name>
    <filename>structtatami__chunked_1_1CustomSparseChunkedOptions.html</filename>
    <member kind="variable">
      <type>size_t</type>
      <name>maximum_cache_size</name>
      <anchorfile>structtatami__chunked_1_1CustomSparseChunkedOptions.html</anchorfile>
      <anchor>a7590657e09265b2af90f3b493c82fc6e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>require_minimum_cache</name>
      <anchorfile>structtatami__chunked_1_1CustomSparseChunkedOptions.html</anchorfile>
      <anchor>a35d3fd4f3063bfb2f8a5808c8e41b8f7</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami_chunked::LruSlabCache</name>
    <filename>classtatami__chunked_1_1LruSlabCache.html</filename>
    <templarg>typename Id_</templarg>
    <templarg>class Slab_</templarg>
    <member kind="function">
      <type></type>
      <name>LruSlabCache</name>
      <anchorfile>classtatami__chunked_1_1LruSlabCache.html</anchorfile>
      <anchor>a022fccffb0fcfbe9e399366eb6f22bb1</anchor>
      <arglist>(size_t m=1)</arglist>
    </member>
    <member kind="function">
      <type>const Slab_ &amp;</type>
      <name>find</name>
      <anchorfile>classtatami__chunked_1_1LruSlabCache.html</anchorfile>
      <anchor>a022825f5749bba03b112fbe72d98ddd5</anchor>
      <arglist>(Id_ id, Cfunction_ create, Pfunction_ populate)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami_chunked::MockSimpleDenseChunk</name>
    <filename>structtatami__chunked_1_1MockSimpleDenseChunk.html</filename>
    <member kind="typedef">
      <type>double</type>
      <name>value_type</name>
      <anchorfile>structtatami__chunked_1_1MockSimpleDenseChunk.html</anchorfile>
      <anchor>a88f07a7b34a0d7a0f9fbd296467f24aa</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>std::vector&lt; value_type &gt;</type>
      <name>Workspace</name>
      <anchorfile>structtatami__chunked_1_1MockSimpleDenseChunk.html</anchorfile>
      <anchor>a0a377b07a002320195e9a0be23c05e97</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1MockSimpleDenseChunk.html</anchorfile>
      <anchor>a3302828143effaa3b331331c64565d3a</anchor>
      <arglist>(Index_ secondary_start, Index_ secondary_length, Workspace &amp;work, value_type *output, size_t stride) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1MockSimpleDenseChunk.html</anchorfile>
      <anchor>a22b1030cdf9149c775b7a45677a73086</anchor>
      <arglist>(const std::vector&lt; Index_ &gt; &amp;secondary_indices, Workspace &amp;work, value_type *output, size_t stride) const</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static constexpr bool</type>
      <name>use_subset</name>
      <anchorfile>structtatami__chunked_1_1MockSimpleDenseChunk.html</anchorfile>
      <anchor>ab8b6933c5a77ea1d9e29ee4d91413b3d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami_chunked::MockSimpleSparseChunk</name>
    <filename>structtatami__chunked_1_1MockSimpleSparseChunk.html</filename>
    <class kind="struct">tatami_chunked::MockSimpleSparseChunk::Workspace</class>
    <member kind="typedef">
      <type>double</type>
      <name>value_type</name>
      <anchorfile>structtatami__chunked_1_1MockSimpleSparseChunk.html</anchorfile>
      <anchor>ad6c1f27da026f8d947159884d5e4cd1f</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>int</type>
      <name>index_type</name>
      <anchorfile>structtatami__chunked_1_1MockSimpleSparseChunk.html</anchorfile>
      <anchor>a8aac93a4ea386e10341de5a0626fc6a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1MockSimpleSparseChunk.html</anchorfile>
      <anchor>a0eda6522f83a7e849148597ed8e3d150</anchor>
      <arglist>(Index_ secondary_start, Index_ secondary_length, Workspace &amp;work, std::vector&lt; std::vector&lt; value_type &gt; &gt; &amp;output_values, std::vector&lt; std::vector&lt; index_type &gt; &gt; &amp;output_indices, index_type shift) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1MockSimpleSparseChunk.html</anchorfile>
      <anchor>a91cf20517832532e1616769a14615cf1</anchor>
      <arglist>(const std::vector&lt; Index_ &gt; &amp;secondary_indices, Workspace &amp;work, std::vector&lt; std::vector&lt; value_type &gt; &gt; &amp;output_values, std::vector&lt; std::vector&lt; index_type &gt; &gt; &amp;output_indices, index_type shift) const</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static constexpr bool</type>
      <name>use_subset</name>
      <anchorfile>structtatami__chunked_1_1MockSimpleSparseChunk.html</anchorfile>
      <anchor>aacc656011869b48792e99e90b393110a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami_chunked::MockSubsetDenseChunk</name>
    <filename>structtatami__chunked_1_1MockSubsetDenseChunk.html</filename>
    <member kind="typedef">
      <type>double</type>
      <name>value_type</name>
      <anchorfile>structtatami__chunked_1_1MockSubsetDenseChunk.html</anchorfile>
      <anchor>af44118c0e432a492297ad1f9cea0ac4d</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>std::vector&lt; value_type &gt;</type>
      <name>Workspace</name>
      <anchorfile>structtatami__chunked_1_1MockSubsetDenseChunk.html</anchorfile>
      <anchor>a81fdb711908f7ab2a0d22306dc369b2d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1MockSubsetDenseChunk.html</anchorfile>
      <anchor>a0b918bb9dfd60289ae743cbd61efba7b</anchor>
      <arglist>(Index_ primary_start, Index_ primary_length, Index_ secondary_start, Index_ secondary_length, Workspace &amp;work, value_type *output, size_t stride) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1MockSubsetDenseChunk.html</anchorfile>
      <anchor>a5ccdaa43b292223a7acdb8eaa24d3603</anchor>
      <arglist>(Index_ primary_start, Index_ primary_length, const std::vector&lt; Index_ &gt; &amp;secondary_indices, Workspace &amp;work, value_type *output, size_t stride) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1MockSubsetDenseChunk.html</anchorfile>
      <anchor>abbe1d50ad042df17cfdd30ea2a8d0ecd</anchor>
      <arglist>(const std::vector&lt; Index_ &gt; &amp;primary_indices, Index_ secondary_start, Index_ secondary_length, Workspace &amp;work, value_type *output, size_t stride) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1MockSubsetDenseChunk.html</anchorfile>
      <anchor>a21e900985cf3160e816e6eff5b715dfc</anchor>
      <arglist>(const std::vector&lt; Index_ &gt; &amp;primary_indices, const std::vector&lt; Index_ &gt; &amp;secondary_indices, Workspace &amp;work, value_type *output, size_t stride) const</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static constexpr bool</type>
      <name>use_subset</name>
      <anchorfile>structtatami__chunked_1_1MockSubsetDenseChunk.html</anchorfile>
      <anchor>a79e60970d114d27efa8db67042b89b4c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami_chunked::MockSubsetSparseChunk</name>
    <filename>structtatami__chunked_1_1MockSubsetSparseChunk.html</filename>
    <class kind="struct">tatami_chunked::MockSubsetSparseChunk::Workspace</class>
    <member kind="typedef">
      <type>double</type>
      <name>value_type</name>
      <anchorfile>structtatami__chunked_1_1MockSubsetSparseChunk.html</anchorfile>
      <anchor>ab4b8726956944e3d23a847873467733d</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>int</type>
      <name>index_type</name>
      <anchorfile>structtatami__chunked_1_1MockSubsetSparseChunk.html</anchorfile>
      <anchor>ae57b40750c649c924a5c8b7d5b564c98</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1MockSubsetSparseChunk.html</anchorfile>
      <anchor>ad8d849ed25414aad22e85df501ee9020</anchor>
      <arglist>(Index_ primary_start, Index_ primary_length, Index_ secondary_start, Index_ secondary_length, Workspace &amp;work, std::vector&lt; std::vector&lt; value_type &gt; &gt; &amp;output_values, std::vector&lt; std::vector&lt; index_type &gt; &gt; &amp;output_indices, index_type shift) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1MockSubsetSparseChunk.html</anchorfile>
      <anchor>add2224603efe08526d6ccc930960a922</anchor>
      <arglist>(Index_ primary_start, Index_ primary_length, const std::vector&lt; Index_ &gt; &amp;secondary_indices, Workspace &amp;work, std::vector&lt; std::vector&lt; value_type &gt; &gt; &amp;output_values, std::vector&lt; std::vector&lt; index_type &gt; &gt; &amp;output_indices, index_type shift) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1MockSubsetSparseChunk.html</anchorfile>
      <anchor>ab2a62c44d054c8c66fa6219079340ec4</anchor>
      <arglist>(const std::vector&lt; Index_ &gt; &amp;primary_indices, Index_ secondary_start, Index_ secondary_length, Workspace &amp;work, std::vector&lt; std::vector&lt; value_type &gt; &gt; &amp;output_values, std::vector&lt; std::vector&lt; index_type &gt; &gt; &amp;output_indices, index_type shift) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1MockSubsetSparseChunk.html</anchorfile>
      <anchor>aeed5a614e70643585ee2bcd1b3769619</anchor>
      <arglist>(const std::vector&lt; Index_ &gt; &amp;primary_indices, const std::vector&lt; Index_ &gt; &amp;secondary_indices, Workspace &amp;work, std::vector&lt; std::vector&lt; value_type &gt; &gt; &amp;output_values, std::vector&lt; std::vector&lt; index_type &gt; &gt; &amp;output_indices, index_type shift) const</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static constexpr bool</type>
      <name>use_subset</name>
      <anchorfile>structtatami__chunked_1_1MockSubsetSparseChunk.html</anchorfile>
      <anchor>ae184bc71b5796bc03d981cc74bad7a26</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami_chunked::OracleSlabCache</name>
    <filename>classtatami__chunked_1_1OracleSlabCache.html</filename>
    <templarg>typename Id_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>class Slab_</templarg>
    <member kind="function">
      <type></type>
      <name>OracleSlabCache</name>
      <anchorfile>classtatami__chunked_1_1OracleSlabCache.html</anchorfile>
      <anchor>a920b4c5979177df0df54c88cc77f399a</anchor>
      <arglist>(std::shared_ptr&lt; const tatami::Oracle&lt; Index_ &gt; &gt; ora, size_t num_slabs)</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>next</name>
      <anchorfile>classtatami__chunked_1_1OracleSlabCache.html</anchorfile>
      <anchor>a368d924ebc4d47daf8f0766d9c0a2438</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; const Slab_ *, Index_ &gt;</type>
      <name>next</name>
      <anchorfile>classtatami__chunked_1_1OracleSlabCache.html</anchorfile>
      <anchor>a8a7fd554334053ebef56bdcb5d6c0321</anchor>
      <arglist>(Ifunction_ identify, Cfunction_ create, Pfunction_ populate)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami_chunked::SimpleDenseChunkWrapper</name>
    <filename>structtatami__chunked_1_1SimpleDenseChunkWrapper.html</filename>
    <templarg>class Blob_</templarg>
  </compound>
  <compound kind="struct">
    <name>tatami_chunked::SimpleSparseChunkWrapper</name>
    <filename>structtatami__chunked_1_1SimpleSparseChunkWrapper.html</filename>
    <templarg>class Blob_</templarg>
  </compound>
  <compound kind="class">
    <name>tatami_chunked::SubsettedOracleSlabCache</name>
    <filename>classtatami__chunked_1_1SubsettedOracleSlabCache.html</filename>
    <templarg>typename Id_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>class Slab_</templarg>
    <member kind="function">
      <type></type>
      <name>SubsettedOracleSlabCache</name>
      <anchorfile>classtatami__chunked_1_1SubsettedOracleSlabCache.html</anchorfile>
      <anchor>ac642e224ec611e618a9c2a8fa9b6c4bc</anchor>
      <arglist>(std::shared_ptr&lt; const tatami::Oracle&lt; Index_ &gt; &gt; ora, size_t num_slabs)</arglist>
    </member>
    <member kind="function">
      <type>Index_</type>
      <name>next</name>
      <anchorfile>classtatami__chunked_1_1SubsettedOracleSlabCache.html</anchorfile>
      <anchor>ab861a229c86fab7212a18b05a01fb9a4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; const Slab_ *, Index_ &gt;</type>
      <name>next</name>
      <anchorfile>classtatami__chunked_1_1SubsettedOracleSlabCache.html</anchorfile>
      <anchor>a835433da5c6e146e76f2e50a43d6ce28</anchor>
      <arglist>(Ifunction_ identify, Cfunction_ create, Pfunction_ populate)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami_chunked::SubsettedOracleSlabSubset</name>
    <filename>structtatami__chunked_1_1SubsettedOracleSlabSubset.html</filename>
    <templarg>typename Index_</templarg>
    <member kind="variable">
      <type>SubsettedOracleSelection</type>
      <name>selection</name>
      <anchorfile>structtatami__chunked_1_1SubsettedOracleSlabSubset.html</anchorfile>
      <anchor>a5ef9c05815bf43353e56597543280111</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Index_</type>
      <name>block_start</name>
      <anchorfile>structtatami__chunked_1_1SubsettedOracleSlabSubset.html</anchorfile>
      <anchor>a9bce81b2fe9e82ffeabed184d757bf52</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Index_</type>
      <name>block_length</name>
      <anchorfile>structtatami__chunked_1_1SubsettedOracleSlabSubset.html</anchorfile>
      <anchor>a8f630419cd5d16f450dc7dc2f7d021f1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Index_ &gt;</type>
      <name>indices</name>
      <anchorfile>structtatami__chunked_1_1SubsettedOracleSlabSubset.html</anchorfile>
      <anchor>ac5526368f24a84bf4240b3b50bef1972</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::unordered_map&lt; Index_, Index_ &gt;</type>
      <name>mapping</name>
      <anchorfile>structtatami__chunked_1_1SubsettedOracleSlabSubset.html</anchorfile>
      <anchor>a73221964703115aa90ee69608ed133f6</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami_chunked::TypicalSlabCacheWorkspace</name>
    <filename>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</filename>
    <templarg>bool oracle_</templarg>
    <templarg>bool subset_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>class Slab_</templarg>
    <member kind="typedef">
      <type>std::conditional&lt; oracle_, typenamestd::conditional&lt; subset_, SubsettedOracleSlabCache&lt; Index_, Index_, Slab_ &gt;, OracleSlabCache&lt; Index_, Index_, Slab_ &gt; &gt;::type, LruSlabCache&lt; Index_, Slab_ &gt; &gt;::type</type>
      <name>SlabCache</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>a6197681c15be36666c0f84d59cbbaf6d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TypicalSlabCacheWorkspace</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>ac0a1044ead9aa4efeee929a554bebf1f</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TypicalSlabCacheWorkspace</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>a37ee1eceb475726fa3164040d3bdf036</anchor>
      <arglist>(Index_ primary_length, Index_ secondary_length, size_t cache_size_in_elements, bool require_minimum_cache, tatami::MaybeOracle&lt; oracle_, Index_ &gt; oracle)</arglist>
    </member>
    <member kind="variable">
      <type>size_t</type>
      <name>slab_size_in_elements</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>a9dbdcf92bd2886c49cc17c3d1c2f4199</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>size_t</type>
      <name>num_slabs_in_cache</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>a0a3b4426fcfe4d5b5bef07e4dd83b86a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>SlabCache</type>
      <name>cache</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>a6a786e1be6531f07dd5f8718ca4526d9</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami_chunked::MockSimpleSparseChunk::Workspace</name>
    <filename>structtatami__chunked_1_1MockSimpleSparseChunk_1_1Workspace.html</filename>
  </compound>
  <compound kind="struct">
    <name>tatami_chunked::MockSubsetSparseChunk::Workspace</name>
    <filename>structtatami__chunked_1_1MockSubsetSparseChunk_1_1Workspace.html</filename>
  </compound>
  <compound kind="namespace">
    <name>tatami_chunked</name>
    <filename>namespacetatami__chunked.html</filename>
    <class kind="class">tatami_chunked::CustomDenseChunkedMatrix</class>
    <class kind="struct">tatami_chunked::CustomDenseChunkedOptions</class>
    <class kind="class">tatami_chunked::CustomSparseChunkedMatrix</class>
    <class kind="struct">tatami_chunked::CustomSparseChunkedOptions</class>
    <class kind="class">tatami_chunked::LruSlabCache</class>
    <class kind="struct">tatami_chunked::MockSimpleDenseChunk</class>
    <class kind="struct">tatami_chunked::MockSimpleSparseChunk</class>
    <class kind="struct">tatami_chunked::MockSubsetDenseChunk</class>
    <class kind="struct">tatami_chunked::MockSubsetSparseChunk</class>
    <class kind="class">tatami_chunked::OracleSlabCache</class>
    <class kind="struct">tatami_chunked::SimpleDenseChunkWrapper</class>
    <class kind="struct">tatami_chunked::SimpleSparseChunkWrapper</class>
    <class kind="class">tatami_chunked::SubsettedOracleSlabCache</class>
    <class kind="struct">tatami_chunked::SubsettedOracleSlabSubset</class>
    <class kind="struct">tatami_chunked::TypicalSlabCacheWorkspace</class>
    <member kind="enumeration">
      <type></type>
      <name>SubsettedOracleSelection</name>
      <anchorfile>namespacetatami__chunked.html</anchorfile>
      <anchor>ad8a9a6497dbbe538be338a3b0fc3bcea</anchor>
      <arglist></arglist>
      <enumvalue file="namespacetatami__chunked.html" anchor="ad8a9a6497dbbe538be338a3b0fc3bceaaba7de5bc6888294e5884b024a4c894f1">FULL</enumvalue>
      <enumvalue file="namespacetatami__chunked.html" anchor="ad8a9a6497dbbe538be338a3b0fc3bceaa4d34f53389ed7f28ca91fc31ea360a66">BLOCK</enumvalue>
      <enumvalue file="namespacetatami__chunked.html" anchor="ad8a9a6497dbbe538be338a3b0fc3bceaacb4ae3b37047fb4b2c0d16f8bf84f076">INDEX</enumvalue>
    </member>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>tatami bindings for chunked matrices</title>
    <filename>index.html</filename>
    <docanchor file="index.html">md__2github_2workspace_2README</docanchor>
  </compound>
</tagfile>
